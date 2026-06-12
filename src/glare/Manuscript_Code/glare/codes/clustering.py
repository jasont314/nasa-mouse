import numpy as np
from tqdm import tqdm
# Clustering
from sklearn.mixture import GaussianMixture
from sklearn.cluster import AgglomerativeClustering
import hdbscan
from sklearn.cluster import SpectralClustering


class GLARECluster:
    def __init__(self, df, representation, location):
        self.df = df
        self.representation = representation
        self.location = location

    def gmm_cluster(self):
        # numer of cluster change depending on a location
        num_cluster = 20 if self.location == 'FLT' else 25
        # GMM # value for max_iter can be changed after tuning
        gmm = GaussianMixture(n_components=num_cluster, random_state=2024, n_init=num_cluster,
                              max_iter=100).fit_predict(self.representation)
        self.df['gmm'] = gmm
        # # Change the dtype of the cluster label
        # df['gmm'] = df['gmm'].apply(str)

        return self.df

    def hdbscan_cluster(self):
        # HDBSCAN doesn't require number of clusters # Set hyperparameter # Can be changed after tuning
        hp_cs = 60 if self.location == 'FLT' else 50
        hdbscan_result = hdbscan.HDBSCAN(min_cluster_size=hp_cs, min_samples= int(hp_cs / 2),
                                         cluster_selection_method='leaf',
                                         prediction_data=True).fit(self.representation)
        # membership vector for soft clustering
        soft_clusters_base = hdbscan.all_points_membership_vectors(hdbscan_result)
        # Assign the highest probability
        soft_c = [np.argmax(x) for x in soft_clusters_base]
        self.df['hdbscan'] = soft_c
        # # Change the dtype of the cluster label
        # df['hdbscan'] = df['hdbscan'].apply(str)

        return self.df

    def spectral_cluster(self):
        # numer of cluster change depending on a location
        num_cluster = 25 if self.location == 'FLT' else 20
        # Spectral Clustering
        sc = SpectralClustering(n_clusters=num_cluster, affinity='nearest_neighbors', n_jobs=-1,
                                random_state=2024).fit(self.representation)
        self.df['spectral'] = sc.labels_
        # # Change the dtype of the cluster label
        # df['spectral'] = df['spectral'].apply(str)

        return self.df

    def eac(self):
        """
            Combine multiple clustering labels using Evidence Accumulation Clustering (EAC).

            Args:
                label_df (df): A pandas dataframe, containing cluter labels from gmm, HDBSCAN, and spectral

            Returns:
                label_df: Dataframe with the consensus clustering obtained by EAC.
        """
        # Run base clustering algorithms to get all labels
        with tqdm(total=3, desc="Running Base Clustering...") as pbar:
            # Run GMM
            self.df = self.gmm_cluster()
            pbar.update(1)
            # Run HDBSCAN
            self.df = self.hdbscan_cluster()
            pbar.update(1)
            # Run Spectral clustering
            labeled_df = self.spectral_cluster()
            pbar.update(1)
        # If using other clustering algorithm, change the column name accordingly
        cluster_arr = [np.array(labeled_df['gmm']), np.array(labeled_df['hdbscan']),
                       np.array(labeled_df['spectral'])]
        # number of object / loci
        n_objects = len(cluster_arr[0])

        # Initialize similarity matrix
        similarity_matrix = np.zeros((n_objects, n_objects))

        # Accumulate evidence across clusterings by consensus voting
        for clustering in tqdm(cluster_arr, desc='Starting Ensemble Clustering...'):
            for i in range(n_objects):
                for j in range(i + 1, n_objects):
                    if clustering[i] == clustering[j]:
                        similarity_matrix[i, j] += 1
                        similarity_matrix[j, i] += 1

        # Obtain consensus clustering using hierarchical clustering
        consensus_clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=0,
                                                       linkage='average', affinity='precomputed')
        consensus_clustering.fit(1 - similarity_matrix)
        # Assign consensus clustering
        labeled_df['consensus'] = consensus_clustering.labels_

        return labeled_df
