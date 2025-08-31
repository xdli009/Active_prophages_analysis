import pandas as pd
import os

def process_blast_results(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None)
    
    grouped = df.groupby(1)
    
    result_dict = {}

    for group_name, group_df in grouped:
        filtered_df = group_df[group_df[3] >= 5000].copy()
        filtered_df[8], filtered_df[9] = filtered_df[[8, 9]].min(axis=1), filtered_df[[8, 9]].max(axis=1)
        sorted_df = filtered_df.sort_values(by=[8, 9])
        ranges = []
        for idx, row in sorted_df.iterrows():
            start, end = row[8], row[9]
            found = False
            for r in ranges:
                if abs(r[0] - start) <= 1000 and abs(r[1] - end) <= 1000:
                    r[2] += 1
                    query_range = row[0] + ":" + str(row[6]) + "-" + str(row[7])
                    r[3].append(query_range)
                    found = True
                    break
            if not found:
                query_range = row[0] + ":" + str(row[6]) + "-" + str(row[7])
                ranges.append([start, end, 1, [query_range]])
        
        result_dict[group_name] = ranges
    
    return result_dict

dict_length = {}
fr_length = open("/beegfs/home/lxd/6_prophage/4.prophage_genetics/1.classification/0.normal_phage_info/genome_length/phage_normal_length.txt","r").read().splitlines()
for line in fr_length:
    phageid = line.split("\t")[0]
    length = line.split("\t")[1]
    dict_length[phageid] = length
file_path = 'all_to_positive.m8'
result_dict = process_blast_results(file_path)
for group, ranges in result_dict.items():
    mkdir = f"mkdir ./results/{group}"
    os.system(mkdir)
    path = f"./results/{group}/{group}_ranges.txt"
    fw_range = open(path,"w")
    length = dict_length[group]
    for r in ranges:
        phageids = ",".join(r[3])
        fw_range.write(group + "\t" + length + "\t" + str(r[0]) + "\t" + str(r[1]) + "\t" + str(r[2]) + "\t" + phageids + "\n")

