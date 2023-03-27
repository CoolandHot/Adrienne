import pandas as pd
import glob

preceeding_path = "diff_expr/output"

for filename in glob.glob(f"{preceeding_path}/*.csv"):
    df = pd.read_csv(filename).dropna()
    df.rename(columns={df.columns[0]: 'gene'}, inplace=True)
    df.sort_values(by='log2FoldChange', ascending=False, inplace=True)
    df[df['padj'] < 0.05].to_csv(filename, index=False)
