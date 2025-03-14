import pandas as pd
import os 

files_dict = {}
for file in snakemake.input.files:
    sample = file.split("/")[-1].split("_filtered")[0]
    files_dict[sample] = file

# Read TSV file
tsv_df = pd.read_csv(snakemake.input.setup, sep="\t",index_col=None)
group = snakemake.wildcards.group
group_df = tsv_df[tsv_df['group'] == group]

sample = group_df['id']
tissue = group_df['tissue']
df = pd.DataFrame(columns=["file_name","cap_flag","merge_priority(start,junctions,end)","source_name"])
for i in range(0,len(sample)):
    df.loc[i] = [files_dict[sample.iloc[i]],"capped","1,1,1",tissue.iloc[i]]
# Save the dataframe to a file
df.to_csv(snakemake.output[0],sep="\t",index=False,header=False)

