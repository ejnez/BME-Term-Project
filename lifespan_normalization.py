import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression

def estimate_and_normalize_lifespan(df, output_file_path):
    """
    Estimates missing Average longevity using Maximum longevity and normalizes lifespan.
    Normalization is based on a regression-derived factor (alpha) that adjusts Average longevity.
    
    :param df: DataFrame with 'Maximum longevity (yrs)' and 'Average longevity' columns.
    :param output_file_path: Path for the new CSV file to save the results.
    :return: Updated DataFrame.
    """
    train_data = df.dropna(subset=["Average longevity"])

    X_train = train_data[["Maximum longevity (yrs)"]].values 
    y_train = train_data["Average longevity"].values 

    model = LinearRegression()
    model.fit(X_train, y_train)

    alpha = model.coef_[0]  

    missing_avg = df["Average longevity"].isna()
    df.loc[missing_avg, "Average longevity"] = model.predict(df.loc[missing_avg, ["Maximum longevity (yrs)"]])

    df["Normalized lifespan"] = df["Average longevity"] + \
                                (df["Maximum longevity (yrs)"] - df["Average longevity"]) * alpha

    result_df = df[["Organism Name", "Maximum longevity (yrs)", "Average longevity", "Normalized lifespan"]]

    result_df.to_csv(output_file_path, index=False)

    return result_df

file_path = "large_res_data.csv"  
output_file_path = "large_res_data_normalized.csv"  

df = pd.read_csv(file_path)  

df = estimate_and_normalize_lifespan(df, output_file_path)

