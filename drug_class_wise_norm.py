import argparse
import pandas as pd

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--norm", type = str, required = True, help = "rpoB normalization file of ARG (required)")
    parser.add_argument("-o", "--out", type = str, help = "User provided output file name (required)")

    args = parser.parse_args()

    original_rpoB_norm = args.norm
    out_fname = args.out

    df = pd.read_csv(original_rpoB_norm, sep='\t')

    grouped_df = df.groupby('drug')['rpob_Normalization'].sum().reset_index()

    output_file = out_fname + "_drug_wise_rpoB_norm.tsv"
    # Write the DataFrame to a TSV file
    grouped_df.to_csv(output_file, sep='\t', index=False)
