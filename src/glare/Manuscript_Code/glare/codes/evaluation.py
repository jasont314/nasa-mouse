import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, accuracy_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import  KFold
from sklearn.manifold import trustworthiness


def trustworthiness_score(source_df, representations):
    print(round(trustworthiness(source_df.iloc[:,1:-1], representations, n_neighbors=500),3))


def knn_classifier(eval_clusters, representations):
    kf = KFold(n_splits=5, shuffle=True, random_state=2024)
    acc = []
    # enumerate splits
    for train, test in kf.split(eval_clusters):
        X, y = representations, eval_clusters
        X_train, X_test, y_train, y_test = X[train], X[test], y[train], y[test]
        # KNN accuracy # Select appropriate hyperparameter depending on your dataset
        knn = KNeighborsClassifier(n_neighbors=500, metric='cosine', n_jobs=-1)
        # Fit the model on the training data
        knn.fit(X_train, y_train)
        # Make predictions on the test data
        y_pred = knn.predict(X_test)
        acc.append(accuracy_score(y_test, y_pred))
    # Print the result of 5-fold CV of KNN classifier
    print(round(np.mean(acc), 3))


def kmeans_silhouette(representations):
    # k-means for evaluation # Same number of cluster for all representation
    # number set based on with GO analysis on eac clustering result
    eval_kmeans = KMeans(n_clusters=15, random_state=2024)
    eval_clusters = eval_kmeans.fit_predict(representations)
    silhouette_sc = silhouette_score(representations, eval_clusters)
    print(f"Silhouette Score: {silhouette_sc:.4f}")

    return eval_clusters


if __name__ == "__main__":
    # Read all representation df
    # FLT
    source_flt = pd.read_csv('GLDS120_restructured_FLT.csv')
    pca_flt = pd.read_csv('pca_df_flt.csv') # Make it numpy # Use `np.array(pca_df[['x','y']])`
    tsne_flt = pd.read_csv('tsne_df_flt.csv') # Make it numpy # Use `np.array(tsne_df[['x','y']])`
    umap_flt = pd.read_csv('umap_df_flt.csv') # Make it numpy # Use `np.array(umap_df[['x','y']])`
    sae_flt =np.load('SAE_FLT_represntation.npy')
    FTsae_flt = np.load('FTSAE_FLT_representation.npy')
    # Wrap dfs
    FLT_df_lst = [np.array(pca_flt[['x','y']]), np.array(tsne_flt[['x','y']]), np.array(umap_flt[['x','y']]), sae_flt, FTsae_flt]
    # Evaluation run 
    for df in FLT_df_lst:
        clusters = kmeans_silhouette(df) # Use different dataset for their results
        knn_classifier(clusters, df)
        trustworthiness(source_flt, df)
    # GC
    source_gc = pd.read_csv('GLDS120_restructured_GC.csv')
    pca_gc = pd.read_csv('pca_df_gc.csv') # Make it numpy # Use `np.array(pca_df[['x','y']])`
    tsne_gc = pd.read_csv('tsne_df_gc.csv') # Make it numpy # Use `np.array(tsne_df[['x','y']])`
    umap_gc = pd.read_csv('umap_df_gc.csv') # Make it numpy # Use `np.array(umap_df[['x','y']])`
    sae_gc =np.load('SAE_GC_represntation.npy')
    FTsae_gc = np.load('FTSAE_GC_representation.npy')
    # Wrap dfs
    GC_df_lst = [np.array(pca_gc[['x','y']]), np.array(tsne_gc[['x','y']]), np.array(umap_gc[['x','y']]), sae_gc, FTsae_gc]
    # Evaluation run 
    for df in GC_df_lst:
        clusters = kmeans_silhouette(df) # Use different dataset for their results
        knn_classifier(clusters, df)
        trustworthiness(source_gc, df)






