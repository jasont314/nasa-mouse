import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import matplotlib
from scipy.stats import zscore
from sklearn.neighbors import LocalOutlierFactor
import re

matplotlib.use('Agg')  # Do not print the plot
# Initialize kmeans parameters
kmeans_kwargs = {
    "init": "k-means++",
    "n_init": 10,
    "random_state": 1,
}

def restructure_data(
    df: pd.DataFrame,
    condition_keywords: list,
    first_keywords: list,
    secondary_keywords: list = None,
    replicate_identifier: str = "Rep",
    id_col: str = "gene_id"
) -> pd.DataFrame:
    """
    Generalized restructuring of GLDS data using data melting approach.

    Args:
        df (pd.DataFrame): Input wide-format DataFrame.
        condition_keywords (list): Main experimental conditions (e.g. ['Flight', 'Ground control']).
        first_keywords (list): Any other factors than FLT/GC suc as Genotypes (e.g. ['Col-0', 'WS', 'elp2-5']).
        secondary_keywords (list): Additional factors that was used in the experiment (e.g. ['Leaves', 'Roots']).
        replicate_identifier (str): Pattern to extract replicate info
        id_col (str): Name of gene ID column.

    Returns:
        pd.DataFrame: Restructured data via data melting.
    """
    expr_cols = [col for col in df.columns if col != id_col]
    metadata = []

    for col in expr_cols:
        # Flatten column if it's a list 
        col_str = col if isinstance(col, str) else col[0]
        components = re.split(r'[,_]', col_str)  # Split by comma or underscore # Most case in GLDS/OSD datasets

        cond = next((c for c in condition_keywords if any(c.lower() in comp.lower() for comp in components)), None)
        genotype = next((g for g in first_keywords if any(g.lower() in comp.lower() for comp in components)), None)

        secondary = None
        if secondary_keywords:
            secondary = next((s for s in secondary_keywords if any(s.lower() in comp.lower() for comp in components)), None)

        rep = None
        for comp in components:
            if replicate_identifier.lower() in comp.lower():
                match = re.search(r'(\d+)', comp)
                if match:
                    rep = f"{replicate_identifier}{match.group(1)}"
                    break

        if not all([cond, genotype, rep]):
            raise ValueError(f"Could not parse all metadata from column: '{col_str}'")

        label_parts = [rep, genotype]
        if secondary:
            label_parts.append(secondary)
        label = "_".join(label_parts)

        metadata.append({
            "original_col": col,
            "condition": cond,
            "label": label
        })

    meta_df = pd.DataFrame(metadata)

    # Melt data
    long_df = df.melt(
        id_vars=[id_col],
        value_vars=meta_df['original_col'].tolist(),
        var_name='original_col',
        value_name='expression'
    )

    # Merge with metadata
    long_df = long_df.merge(meta_df, on='original_col')

    # Pivot to wide
    wide_df = long_df.pivot_table(
        index=[id_col, 'condition'],
        columns='label',
        values='expression'
    ).reset_index()

    # Rename condition to Location (for output consistency)
    wide_df.rename(columns={'condition': 'Location'}, inplace=True)

    # Optional column ordering
    all_labels = sorted([col for col in wide_df.columns if col not in [id_col, 'Location']])
    final_df = wide_df[[id_col] + all_labels + ['Location']]

    return final_df


def concat_df(representation, df, dimension):
    # concat representation coordinates with the nc_df
    if dimension == 2:
        coordinates_df = pd.DataFrame(representation, columns=['x', 'y'])
    elif dimension == 3:
        coordinates_df = pd.DataFrame(representation, columns=['x', 'y', 'z'])
    else:
        raise TypeError("Choose 2 or 3 for dimension")
    coordinates_df['gene_id'] = df['gene_id']
    coordinates_df['Location'] = df['Location']

    return coordinates_df


def initial_pca_elbow(loc_df, location):
    # Get data representation via pca
    # PCA in 3-dimension for 3d-plot
    loc_df_pca = PCA(n_components=3).fit_transform(loc_df.iloc[:, 1:-1])
    # Elbow method for k-means
    # Create list to hold SSE values for each k
    sse_lst = []
    for k in range(1, 11):
        kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
        # You can run it on the raw data without dimensionality reduction but it detects same outliers
        kmeans.fit(loc_df_pca)
        sse_lst.append(kmeans.inertia_)
    # Visualize results
    plt.plot(range(1, 11), sse_lst, markersize=3, marker='o')
    plt.xticks(range(1, 11))
    plt.title("Elbow method")
    plt.xlabel("Number of Clusters")
    plt.ylabel("WCSS")
    plt.savefig('elbow_method_plot_' + location + '.png')

    return loc_df_pca

def initial_kmeans(loc_df, loc_df_pca, location):
    # After elbow method we set num_cluster = 5 for FLT and = 4 for GC
    num_cluster = 5 if location == 'FLT' else 4 # Find the right cluster number for your dataset
    # Cluster with tuned number of cluster
    model = KMeans(n_clusters=num_cluster, **kmeans_kwargs)
    # Predict cluster and make df
    y_clusters = model.fit_predict(loc_df_pca)
    loc_df_pca = pd.DataFrame(loc_df_pca, columns=['x', 'y', 'z'])
    loc_df_pca['cluster'] = y_clusters
    loc_df_pca = loc_df_pca.astype({'cluster': 'str'})
    loc_df_pca['gene_id'] = loc_df['gene_id']

    # export .csv for the preprocessing data, explore how the clusters are distributed
    loc_df_pca.to_csv('outlier_detection_df_' + location + '.csv', index=False)
    # Detected outliers for the GLDS-120 dataset
    outlier_gene_id = ["AT3G41768", "ATMG00020", "AT1G07590"]
    # # Optional outlier detection methods
    # outlier_gene_id = lof_od(loc_df)
    # outlier_gene_id = zscore_pca(loc_df_pca, location)
    # Clean df
    clean_df = loc_df.drop(loc_df[loc_df['gene_id'].isin(outlier_gene_id)].index).reset_index(drop=True)

    return clean_df


def tsne4viz(nc_df, representation, dim):
    # t-SNE, 2-dimension # change / tune perplexity as need
    rep_viz = TSNE(n_components=dim, perplexity=75, random_state=1996, n_jobs=-1,
                   learning_rate='auto').fit_transform(representation)
    # Get full pca_df
    rep_viz_df = concat_df(rep_viz, nc_df, dimension=dim)

    return rep_viz_df


def lof_od(loc_df):
    # Choose appropriate the hyperparameter based on your dataset
    clf = LocalOutlierFactor(n_neighbors=100, contamination=1e-4)
    loc_df['lof'] = clf.fit_predict(loc_df.iloc[:,1:-1])
    # Get outliers based on LOF
    return list(loc_df[loc_df['lof']==-1]['gene_id'])


def zscore_pca(loc_df_pca, location):
    # Get z-score for PC1
    loc_df_pca['z_score'] = zscore(loc_df_pca['x'])
    # Analyze the z-score by using the exported csv # Print out the distribution and choose appropriate threshold for outlier
    loc_df_pca.to_csv('z_score_df_' + location + '.csv', index=False)
    # Get outliers based on z-score analysis
    return list(loc_df_pca[loc_df_pca['z_score']>1]['gene_id'])
