import pandas as pd
import sys

files_dict = {}
if sys.argv[1] is not list:
    file = sys.argv[1]
    sample= file.split("/")[-1].split("_filtered")[0]
    files_dict[sample] = file
else:
    for file in sys.argv[1]:
        sample = file.split("/")[-1].split("_filtered")[0]
        files_dict[sample] = file

# Read TSV file
tsv_df = pd.read_csv(sys.argv[2], sep="\t",index_col=None)
group = sys.argv[3]
group_df = tsv_df[tsv_df['group'] == group]

sample = group_df['id']
tissue = group_df['tissue']
df = pd.DataFrame(columns=["file_name","cap_flag","merge_priority(start,junctions,end)","source_name"])
for i in range(0,len(sample)):
    df.loc[i] = [files_dict[sample.iloc[i]],"capped","1,1,1",tissue.iloc[i]]
# Save the dataframe to a file
df.to_csv(sys.argv[4],sep="\t",index=False,header=False)

