import pandas as pd
import sys

files_dict = {}
for file in sys.argv[0]:
    sample = file.split("/")[-1].split("_filtered")[0]
    files_dict[sample] = file

# Read TSV file
print(f"Reading TSV file: {sys.argv[1]}")
tsv_df = pd.read_csv(sys.argv[1], sep="\t",index_col=None)
group = sys.argv[2]
print(f"Creating filelist for group: {group}")
group_df = tsv_df[tsv_df['group'] == group]

sample = group_df['id']
tissue = group_df['tissue']
df = pd.DataFrame(columns=["file_name","cap_flag","merge_priority(start,junctions,end)","source_name"])
for i in range(0,len(sample)):
    df.loc[i] = [files_dict[sample.iloc[i]],"capped","1,1,1",tissue.iloc[i]]
# Save the dataframe to a file
df.to_csv(sys.argv[3],sep="\t",index=False,header=False)

