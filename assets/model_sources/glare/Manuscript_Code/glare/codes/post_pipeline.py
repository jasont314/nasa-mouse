import numpy as np
import pandas as pd
import json
import shap
from sklearn.preprocessing import StandardScaler
import xgboost
from .utils import restructure_data


class SHAP_PP:
    def __init__(self, df):
        self.df = df
        # Set for training supervised learning on labeled dataset
        self.X = df.iloc[:, 1:-1]
        self.y = df['Location'].values
        self.X_averaged = self.X.groupby(self.X.columns.str[-3:], axis=1).mean()

    def classification_model(self, train_X):
        # Normalize
        scaler = StandardScaler()
        train_X = scaler.fit_transform(train_X)
        # Classifier
        model = xgboost.XGBClassifier(n_estimators=500, max_depth=10, n_jobs=-1, random_state=1996).fit(train_X, self.y)

        return model

    def get_shap_values(self): 
        # compute SHAP values for all independent experiment
        explainer_all= shap.Explainer(self.classification_model(self.X), seed=2025) 
        shap_values_all = explainer_all(self.X)
        # compute SHAP values for average experiemnts by location and genotype
        explainer_avg= shap.Explainer(self.classification_model(self.X_averaged), seed=2025) 
        shap_values_avg = explainer_avg(self.X_averaged)

        return shap_values_all, shap_values_avg
    
    def write_json(self, shap_values, type):
        # Convert SHAP values object to a dictionary
        shap_dict = {
            "values": shap_values.values.tolist(),
            "base_values": shap_values.base_values.tolist(),
            "data": shap_values.data.tolist()
        }

        # Output the dictionary to a JSON file
        with open(f"shap_values_{type}.json", "w") as json_file:
            json.dump(shap_dict, json_file, indent=4)

# Usage Example
if __name__ == "__main__":
    # Get arguments
    # Read GeneLab dataset
    glds120_nc = pd.read_csv('GLDS-120_Normalized_Counts.csv')
    glds120_nc.rename(columns={'Unnamed: 0': 'gene_id'}, inplace=True)
    # Restructured data
    glds120_nc_rs = restructure_data(glds120_nc)
    # Encode the label
    glds120_nc_rs.replace({"Location": {"FLT": 1, "GC": 0}}, inplace=True)
    # Post pipeline analysis: SHAP
    SHAP_analysis = SHAP_PP(glds120_nc_rs)
    shap_values_all, shap_values_avg = SHAP_analysis.get_shap_values()
    # export the values
    SHAP_analysis.write_json(shap_values_avg, 'averaged')
    SHAP_analysis.write_json(shap_values_all, 'all')
    # Beeswarm plot using values # Use avg or all
    shap.plots.beeswarm(shap_values_avg, max_display=20) 
    # Scatter plot using values # Use avg or all
    fig = shap.plots.scatter(
    shap_values_avg, ylabel="SHAP value\nhigher means more likely to be classified as FLT"
    )
    # Get shap values for individual loci at FLT or GC # Use avg or all
    shap.plots.bar(shap_values_avg[38966]) # MYC2 FLT
    shap.plots.bar(shap_values_avg[38967]) # MYC2 GC

