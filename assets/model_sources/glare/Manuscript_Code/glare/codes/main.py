import numpy as np
import pandas as pd
from tqdm import tqdm
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader

# Verification study
from .utils import restructure_data
import xgboost
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, f1_score
from sklearn.model_selection import KFold
# import matplotlib.pyplot as plt
# from sklearn.metrics import RocCurveDisplay, auc, confusion_matrix, classification_report

# Preprocessing
from .utils import initial_pca_elbow, initial_kmeans

# Data representations
# from representation_learning import get_pca, get_tsne, get_umap
from .representation_learning import train_SAE, finetune_SAE_sc, Adapter, SparseAutoEncoder
from scipy.io import mmread
from .utils import tsne4viz
# Clustering
from .clustering import GLARECluster

# Set environment
# os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'max_split_size_mb:64'
torch.cuda.empty_cache()
torch.manual_seed(2023)


def verification(nc, labels = ["FLT", "GC"]):
    # Take restructured data vis discretization
    # Encode the label
    new_nc = nc.copy()
    new_nc["Location"] = new_nc["Location"].replace({labels[0]: 1, labels[1]: 0})
    new_nc = new_nc.infer_objects(copy=False)
    # Set for training supervised learning on labeled dataset
    X = new_nc.iloc[:, 1:-1]
    y = new_nc['Location'].values
    # train_acc, = []
    test_acc_lst, f1_score_lst, roc_auc_lst = [], [], []
    # 5-fold cross validation
    kf = KFold(n_splits=5, random_state=2023, shuffle=True)
    for i, (train_index, test_index) in enumerate(kf.split(X)):
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y[train_index], y[test_index]
        # Training xgboost for classification, with tuned hyperparameter
        xgb_clf = xgboost.XGBClassifier(n_estimators=500, max_depth=10, n_jobs=-1, random_state=1996)
        xgb_clf.fit(X_train, y_train)
        # # Get prediction for the training dataset
        # y_pred_tr = xgb_clf.predict(X_train)
        # Get prediction for the test dataset
        y_pred = xgb_clf.predict(X_test)
        # # For visualization
        # # ROC, AUC Plot
        # ax = plt.gca()
        # rfc_disp = RocCurveDisplay.from_estimator(xgb_clf, X_test, y_test, ax=ax, alpha=0.8)
        # plt.plot([0, 1], [0, 1], color='red', lw=1, linestyle='--')
        # plt.show()
        #
        # # Confusion matrix
        # cfm = confusion_matrix(y_test, y_pred)
        # tp = cfm[0][0]
        # fn = cfm[0][1]
        # fp = cfm[1][0]
        # tn = cfm[1][1]
        #
        # print('Confusion Matrix for Test data:\n')
        # print(cfm)
        #
        # # Misclassification rate
        # rate = (fp + fn) / (tp + tn + fp + fn)
        # print('\nMisclassification rate:', rate)

        # Test accuracy
        test_acc = xgb_clf.score(X_test, y_test) * 100
        test_acc_lst.append(test_acc)
        # f1 score and roc-auc score
        f1_score_lst.append(f1_score(y_test, y_pred, average='binary'))
        roc_auc_lst.append(roc_auc_score(y_test, xgb_clf.predict_proba(X_test)[:,1]))
        # # Overall report
        # print(classification_report(y_test, y_pred))

    return new_nc, test_acc_lst, f1_score_lst, roc_auc_lst


def preprocessing(df):
    # Divide FLT and GC
    FLT_df = df[df['Location'] == 1].reset_index(drop=True)
    GC_df = df[df['Location'] == 0].reset_index(drop=True)
    # Preprocessing
    # Get pca and print elbow method plot for k-menas
    pca_FLT_df = initial_pca_elbow(FLT_df, 'FLT')
    pca_GC_df = initial_pca_elbow(GC_df, 'GC')
    # Perform k-means clustering, and perform intial investigation on exported df
    # Use appropriate outlier detection appropriate to your dataset
    clean_FLT_df = initial_kmeans(FLT_df, pca_FLT_df, 'FLT')
    clean_GC_df = initial_kmeans(GC_df, pca_GC_df, 'GC')

    return clean_FLT_df, clean_GC_df


def SAE_inference(nc_df, device, location):
    # Train SAE
    SAE_model = train_SAE(np.array(nc_df), device, exp_type=location)
    SAE_model.eval()
    # Create a StandardScaler instance
    scaler = StandardScaler()
    # Fit the scaler to your data and transform it
    X = scaler.fit_transform(np.array(nc_df))
    # To tensor
    X = torch.tensor(X, dtype=torch.float32).to(device)
    # Retrieve data representation from bottleneck layer
    with torch.no_grad():
        encoded_data = SAE_model.encoder(X)
    # Convert the encoded data tensor to NumPy array
    SAE_representation = encoded_data.detach().cpu().numpy()

    return SAE_representation


def FTSAE_inference(nc_df, pi_dim, weights, device, ft_epoch):
    # Finetune SAE # Check the log for the total loss while tuning, and pick the best weights accordingly
    FTSAE_model = finetune_SAE_sc(nc_df, pi_dim, weights, device, ft_epoch)
    FTSAE_model.eval()
    # Create a StandardScaler instance
    scaler = StandardScaler()
    # Fit the scaler to your data and transform it
    X = scaler.fit_transform(np.array(nc_df))
    # To tensor
    X = torch.tensor(X, dtype=torch.float32) # .to(device)
    # Use adapter layer
    adapter = Adapter(X.shape[1], pi_dim)
    X = adapter(X).clone().detach()
    # Inference
    input_dim = X.shape[1]
    model = SparseAutoEncoder(input_dim)
    # Initialize the model # Weight from last epoch of fine-tuning step
    finetuned_weights =  f'./weights/finetune_checkpoint_epoch_{ft_epoch}.pth' # './weights/FLT_finetuned.pth' if location == 'FLT' else './weights/GC_finetuned.pth'
    checkpoint = torch.load(finetuned_weights)  # Load the checkpoint file
    # Load the model state dict
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()
    # Retrieve data representation from bottleneck layer
    with torch.no_grad():
        encoded_data = model.encoder(X)

    # Convert the encoded data tensor to NumPy array
    FTSAE_representation = encoded_data.detach().cpu().numpy()

    return FTSAE_representation


if __name__ == "__main__":
    # Set seed
    torch.manual_seed(1996)
    # Set gpu # cuda required, if you are using different gpu change.
    gpu = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    # Read GeneLab dataset
    glds120_nc = pd.read_csv('GLDS-120_Normalized_Counts.csv')
    glds120_nc.rename(columns={'Unnamed: 0': 'gene_id'}, inplace=True)

    # Verification study
    glds120 = restructure_data(
        df=glds120_nc,
        condition_keywords=['FLT', 'GC'],
        first_keywords=['Col-0-PhyD', 'Col-0', 'Ws'],
        secondary_keywords=['Alight', 'dark'],
        replicate_identifier='Rep',
        id_col='gene_id'
    )

    new_glds120, test_acc_lst, f1_score_lst, roc_auc_lst = verification(glds120)
    # Print verification study result
    print('test accuracy:', np.mean(test_acc_lst), '+/-', round(np.std(test_acc_lst), 3))
    print('f1 score:', round(np.mean(f1_score_lst), 5), '+/-', round(np.std(f1_score_lst), 3))
    print('roc-auc:', round(np.mean(roc_auc_lst), 5), '+/-', round(np.std(roc_auc_lst), 3))

    # Preprocessing 
    # While we perform cluter-based outlier detection for preprocessing, choose suitable method for your specific data.
    clean_glds120_FLT, clean_glds120_GC = preprocessing(new_glds120)

    # ## Run other base representation learning algorithm
    # # FLT
    # # PCA, export
    # pca_df_FLT = get_pca(clean_glds120_FLT.iloc[:, 1:-1])
    # pca_df_FLT.to_csv('pca_df_FLT.csv', index=False)
    # # t-SNE, export
    # tsne_df_FLT = get_tsne(clean_glds120_FLT.iloc[:, 1:-1])
    # tsne_df_FLT.to_csv('tsne_df_FLT.csv', index=False)
    # # UMAP, export
    # umap_df_FLT = get_umap(clean_glds120_FLT.iloc[:, 1:-1])
    # umap_df_FLT.to_csv('umap_df_FLT.csv', index=False)
    # # GC
    # # PCA, export
    # pca_df_GC = get_pca(clean_glds120_GC.iloc[:, 1:-1])
    # pca_df_GC.to_csv('pca_df_GC.csv', index=False)
    # # t-SNE, export
    # tsne_df_GC = get_tsne(clean_glds120_GC.iloc[:, 1:-1])
    # tsne_df_GC.to_csv('tsne_df_GC.csv', index=False)
    # # UMAP, export
    # umap_df_GC = get_umap(clean_glds120_GC.iloc[:, 1:-1])
    # umap_df_GC.to_csv('umap_df_GC.csv', index=False)

    # ## Direct SAE training on our target data, normalized counts from GLDS-120, export
    # # FLT
    # SAE_FLT_represntation = SAE_inference(clean_glds120_FLT.iloc[:, 1:-1], gpu, location='FLT')
    # np.save('SAE_FLT_represntation.npy', SAE_FLT_represntation)
    # # GC
    # SAE_GC_represntation = SAE_inference(clean_glds120_GC.iloc[:, 1:-1], gpu, location='GC')
    # np.save('SAE_GC_represntation.npy', SAE_GC_represntation)

    ## Single-Cell pre-training
    # Read Single-cell matrix
    sparse_matrix = mmread('E-CURD-5.aggregated_filtered_normalised_counts.mtx')
    sparse_tensor = torch.sparse.FloatTensor(
        torch.LongTensor([sparse_matrix.row, sparse_matrix.col]),
        torch.FloatTensor(sparse_matrix.data),
        torch.Size(sparse_matrix.shape)
    )
    # Get dimension
    pretrained_input_dim = sparse_tensor.shape[1]  # 3552
    # Pretraining with single cell dataset
    SAE_sc = train_SAE(sparse_tensor.to_dense(), gpu, 'sc_pretrain')

    ## Fine-tuning with GLDS data, export
    # FLT
    FTSAE_FLT_representation = FTSAE_inference(clean_glds120_FLT.iloc[:, 1:-1], pretrained_input_dim,
                                              './weights/sc_shulse_pretrained.pth', gpu, location='FLT')
    np.save('FTSAE_FLT_representation.npy', FTSAE_FLT_representation)
    # GC
    FTSAE_GC_representation = FTSAE_inference(clean_glds120_GC.iloc[:, 1:-1], pretrained_input_dim,
                                              './weights/sc_shulse_pretrained.pth', gpu, location='GC')
    np.save('FTSAE_GC_representation.npy', FTSAE_GC_representation)

    ## For Visaulization
    # FLT
    viz_df_FLT = tsne4viz(clean_glds120_FLT, FTSAE_FLT_representation, dim=2)
    # GC
    viz_df_GC = tsne4viz(clean_glds120_GC, FTSAE_GC_representation, dim=2)

    # Clustering: FLT
    glare_cst_FLT = GLARECluster(viz_df_FLT, FTSAE_FLT_representation, 'FLT')
    # # If you want to run the base clustering separately to get them
    # gmm_df_FLT = glare_cst.gmm_cluster()
    # hdbscan_df_FLT = glare_cst.hdbscan_cluster()
    # spectral_df_FLT = glare_cst.spectral_cluster()
    # Run ensemble clustering
    labeled_df_FLT = glare_cst_FLT.eac()
    # Export
    labeled_df_FLT.to_csv('final_FLT_df.csv', index=False)

    # Clustering: GC
    glare_cst_GC = GLARECluster(viz_df_GC, FTSAE_GC_representation, 'GC')
    # # If you want to run the base clustering separately to get them
    # gmm_df_GC = glare_cst.gmm_cluster()
    # hdbscan_df_GC = glare_cst.hdbscan_cluster()
    # spectral_df_GC = glare_cst.spectral_cluster()
    # Run ensemble clustering
    labeled_df_GC = glare_cst_GC.eac()
    # Export
    labeled_df_GC.to_csv('final_GC_df.csv', index=False)
