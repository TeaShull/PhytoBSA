import os
import pandas as pd

from settings.paths import OUTPUT_DIR

def label_dataframes(chrom, pos, line_name):
    """
    Label dataframes with ground truth data.
    """
    # Iterate through all subdirectories in the output directory
    for subdir, dirs, files in os.walk(OUTPUT_DIR):
        for file in files:
            # Check if the file is a csv that needs to be labeled and contains the line name
            if file.endswith('_all.csv') and line_name in file:
                # Load the dataframe
                df = pd.read_csv(os.path.join(subdir, file))
                # If 'causal' column does not exist, create it with all False
                if 'causal' not in df.columns:
                    df['causal'] = False
                # Label the rows that match the given chromosome and position as True
                df.loc[(df['chrom'] == chrom) & (df['pos'] == pos), 'causal'] = True
                # Save the labeled dataframe
                df.to_csv(os.path.join(subdir, file), index=False)

def analyze_labeled_data():
    """
    Analyze the set of labeled csv tables.
    """
    # Placeholder for analysis results
    analysis_results = {}

    # Iterate through all subdirectories in the output directory
    for subdir, dirs, files in os.walk(OUTPUT_DIR):
        for file in files:
            # Check if the file is a csv that needs to be analyzed
            if file.endswith('_all.csv'):
                # Load the dataframe
                df = pd.read_csv(os.path.join(subdir, file))
                # Analyze the dataframe
                # This is a placeholder for your actual analysis code
                analysis_result = analyze_dataframe(df)
                # Store the analysis result
                analysis_results[file] = analysis_result

    return analysis_results

def analyze_dataframe(df):
    """
    Analyze a dataframe.
    Calculate efficiency, precision and accuracy of the features.
    """
    # Placeholder for analysis results
    analysis_result = {}

    # List of features
    features = ['ratio', 'ratio_yhat', 'G_S', 'G_S_yhat', 'RS_G', 'RS_G_yhat']

    # Calculate efficiency, precision and accuracy for each feature
    for feature in features:
        # Calculate efficiency
        efficiency = np.mean(df[feature + '_percentile'] > 0.95)
        # Calculate precision
        precision = np.mean(df[(df[feature + '_percentile'] > 0.95) & (df['causal'] == True)])
        # Calculate accuracy
        accuracy = np.mean((df[feature + '_percentile'] > 0.95) == df['causal'])
        # Store the results
        analysis_result[feature] = {'efficiency': efficiency, 'precision': precision, 'accuracy': accuracy}

    return analysis_result

