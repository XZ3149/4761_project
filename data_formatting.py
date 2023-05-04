# -*- coding: utf-8 -*-
"""
Created on Sat Apr 22 16:52:58 2023

@author: we609
"""

import csv
import random
import pandas as pd
import math
import numpy as np
import warnings
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)



def read_df(ctcf, h3k9me3, dnase, rna): 
    """
    
    this is the function use to read in files

    """

    
    
    
    ctcf_df = pd.read_csv(ctcf, sep = '\t')
    ctcf_df = ctcf_df[ctcf_df['A549'].notnull()]

    h3k9me3_df = pd.read_csv(h3k9me3, sep = '\t')
    h3k9me3_df = h3k9me3_df[h3k9me3_df['A549'].notnull()]
    
    dnase_df = pd.read_csv(dnase, sep = '\t')
    dnase_df = dnase_df[dnase_df['A549'].notnull()]
    
    rna_df = pd.read_csv(rna, sep = '\t')
    rna_df = rna_df[rna_df['A549'].notnull()] 
    
    
    
    
    ctcf_df = normlization(ctcf_df)
    h3k9me3_df = normlization(h3k9me3_df)
    dnase_df = normlization(dnase_df)
    rna_df = normlization(rna_df)
    
    
    
    return ctcf_df, h3k9me3_df, dnase_df, rna_df





def normlization(df):
    """this is the function use to normalize samples, will use read count and log
    """
    
    
    cell_line = ['A549', 'A673', 'Caco-2', 'GM12878', 'H1', 'H9',
           'HCT116', 'HepG2', 'IMR-90', 'K562', 'MCF-7', 'OCI-LY7', 'PC-3', 'PC-9',
           'Panc1']


    """

    for i in cell_line:
        percentile_99 = np.percentile(df[i].astype('int32'), 75)
        
        df[i] = df[i].astype('int32') /percentile_99
        print( percentile_99)
    """  
        
        
    Sum_dic = {}
    total_sum = 0
    for i in cell_line:
        
        Sum_dic[i] = df[i].sum()
        total_sum += df[i].sum()
        
        
    average_sum = total_sum/len(cell_line)
    
    for i in cell_line:
        
        
        df[i] = np.log(df[i]/(Sum_dic[i]/average_sum) + 1)
        
      
        

    return df









def sample_split(ctcf_df, h3k9me3_df, dnase_df, rna_df):
    """ 
    this function is use to split sample into insample group and out sample group
    """
    
    
    ctcf_df[['Chr', 'Start', 'End', 'PC-3', 'OCI-LY7']].to_csv("outsample_ctcf.bin", index=False)
    h3k9me3_df[['Chr', 'Start', 'End', 'PC-3', 'OCI-LY7']].to_csv("outsample_h3k9me3.bin", index=False)
    dnase_df[['Chr', 'Start', 'End', 'PC-3', 'OCI-LY7']].to_csv("outsample_DNase.bin", index=False)
    rna_df[['Chr', 'Start', 'End', 'PC-3', 'OCI-LY7']].to_csv("outsample_rna.bin", index=False)





    ctcf_df[['Chr', 'Start', 'End', 'A549', 'A673', 'Caco-2', 'GM12878', 'H1', 'H9',
           'HCT116', 'HepG2', 'IMR-90', 'K562', 'MCF-7', 'PC-9',
           'Panc1']].to_csv("insample_ctcf.bin", index=False)


    h3k9me3_df[['Chr', 'Start', 'End', 'A549', 'A673', 'Caco-2', 'GM12878', 'H1', 'H9',
           'HCT116', 'HepG2', 'IMR-90', 'K562', 'MCF-7', 'PC-9',
           'Panc1']].to_csv("insample_h3k9me3.bin", index=False)

    
    dnase_df[['Chr', 'Start', 'End', 'A549', 'A673', 'Caco-2', 'GM12878', 'H1', 'H9',
           'HCT116', 'HepG2', 'IMR-90', 'K562', 'MCF-7', 'PC-9',
           'Panc1']].to_csv("insample_DNase.bin", index=False)


    rna_df[['Chr', 'Start', 'End', 'A549', 'A673', 'Caco-2', 'GM12878', 'H1', 'H9',
           'HCT116', 'HepG2', 'IMR-90', 'K562', 'MCF-7', 'PC-9',
           'Panc1']].to_csv("insample_rna.bin", index=False)





def out_sample_process(ctcf_df, h3k9me3_df, dnase_df, rna_df, samples_size, interval = 500000, random_seed = 123):

    
    samples = dnase_df.columns[3:] # sample name start at line 4;
    num_sample = dnase_df.shape[1] - 3

    random_df =dnase_df.sample(n = int((samples_size/0.9) * 2))
    list_ctcf = []
    list_h3k9me3 = []
    list_rna = []
    list_dnas = []
    list_name = []
    
    count = 0
    for index,row in random_df.iterrows():
        sample_index = random.randrange(0, num_sample - 1)   # skip first three line
        
        
        
        current_sample = samples[sample_index]
        
        temp_dnas = row[current_sample]
        temp_ctcf = information_extration(ctcf_df, current_sample, row['Chr'], row['Start'], interval)
        temp_h3k9me3 = information_extration(h3k9me3_df, current_sample, row['Chr'], row['Start'], interval)
        temp_rna = information_extration(rna_df, current_sample, row['Chr'], row['Start'], interval)
        
        
        if sum(temp_ctcf) >= 0.1 and sum(temp_h3k9me3) >=0.1 and sum(temp_rna) >= 0.1 and\
            len(temp_ctcf) == 1001 and len(temp_h3k9me3) == 1001 and len(temp_rna) == 1001:
        
        
            list_dnas += [[temp_dnas]]
            list_ctcf += [temp_ctcf]
            list_h3k9me3 += [temp_h3k9me3]
            list_rna += [temp_rna]
            
    
            list_name += [f'{row["Chr"]}_{row["Start"]}_{current_sample}']
            count += 1
            if count%500== 0:
                print(f'Train Finished {count}')
    
    
        if count >= ((samples_size/0.9) * 0.1):
            break
    
    
    
    
    
    

    result = pd.DataFrame(list_dnas, index = list_name).transpose()
    result.to_csv("outsample_test_DNase.csv", index=False)
    
    result = pd.DataFrame(list_ctcf, index = list_name).transpose()
    result.to_csv("outsample_test_ctcf.csv", index=False)
    
    
    result = pd.DataFrame(list_h3k9me3, index = list_name).transpose()
    result.to_csv("outsample_test_h3k9me3.csv", index=False)
    
    
    result = pd.DataFrame(list_rna, index = list_name).transpose()
    result.to_csv("outsample_test_rna.csv", index=False)
    






def random_sample_generation(ctcf_df, h3k9me3_df, dnase_df, rna_df, samples_size, interval = 500000, random_seed = 123):

     
    
    samples = dnase_df.columns[3:] # sample name start at line 4;
    num_sample = dnase_df.shape[1] - 3
    
    """
    ctcf_df, h3k9me3_df, dnase_df, rna_df = ctcf_df.sort_values(by = ['Chr', 'Start']), \
        h3k9me3_df.sort_values(by = ['Chr', 'Start']), dnase_df.sort_values(by = ['Chr', 'Start']),\
        rna_df.sort_values(by = ['Chr', 'Start'])
    """
    
    
    
    N = len(dnase_df ) 
    np.random.seed(random_seed) 
    rand_perm = np.random.permutation(N)
    train_idx = rand_perm[:int(np.ceil(0.8 * N))]
    val_idx = rand_perm[int(np.ceil(0.8 * N)):int(np.ceil(0.9 * N))]
    test_idx = rand_perm[int(np.ceil(0.9 * N)):]

    
    
    
    df_train = dnase_df.iloc[train_idx, :].sort_values(by=['Chr', 'Start']).reset_index(drop=True)
    df_val = dnase_df.iloc[val_idx, :].sort_values(by=['Chr', 'Start']).reset_index(drop=True)
    df_test = dnase_df.iloc[test_idx, :].sort_values(by=['Chr', 'Start']).reset_index(drop=True)





    #processing training samples
    random_df =df_train.sample(n = 20 * samples_size)
    list_ctcf = []
    list_h3k9me3 = []
    list_rna = []
    list_dnas = []
    list_name = []
    
    count = 0
    
    for index,row in random_df.iterrows():
        sample_index = random.randrange(0, num_sample - 1)   # skip first three line
        
        
        
        current_sample = samples[sample_index]
        
        temp_dnas = row[current_sample]
        temp_ctcf = information_extration(ctcf_df, current_sample, row['Chr'], row['Start'], interval)
        temp_h3k9me3 = information_extration(h3k9me3_df, current_sample, row['Chr'], row['Start'], interval)
        temp_rna = information_extration(rna_df, current_sample, row['Chr'], row['Start'], interval)
        
        
        
        if sum(temp_ctcf) >= 0.1 and sum(temp_h3k9me3) >=0.1 and sum(temp_rna) >= 0.1 and\
            len(temp_ctcf) == 1001 and len(temp_h3k9me3) == 1001 and len(temp_rna) == 1001:
        
            list_dnas += [[temp_dnas]]
            list_ctcf += [temp_ctcf]
            list_h3k9me3 += [temp_h3k9me3]
            list_rna += [temp_rna]
            
    
            list_name += [f'{row["Chr"]}_{row["Start"]}_{current_sample}']
            count += 1
            if count%500== 0:
                print(f'Train Finished {count}')
    
    
        if count >= samples_size:
            break
    
    
    result = pd.DataFrame(list_dnas, index = list_name).transpose()
    result.to_csv("train_DNase.csv", index=False)
    
    result = pd.DataFrame(list_ctcf, index = list_name).transpose()
    result.to_csv("train_ctcf.csv", index=False)
    
    
    result = pd.DataFrame(list_h3k9me3, index = list_name).transpose()
    result.to_csv("train_h3k9me3.csv", index=False)
    
    
    result = pd.DataFrame(list_rna, index = list_name).transpose()
    result.to_csv("train_rna.csv", index=False)
    
    
    
    
    
    
    
    #processing training samples
    random_df =df_val.sample(n = int((samples_size/0.9) * 2.0))
    list_ctcf = []
    list_h3k9me3 = []
    list_rna = []
    list_dnas = []
    list_name = []
    
    count = 0
    for index,row in random_df.iterrows():
        sample_index = random.randrange(0, num_sample - 1)   # skip first three line
        
        
        
        current_sample = samples[sample_index]
        
        temp_dnas = row[current_sample]
        temp_ctcf = information_extration(ctcf_df, current_sample, row['Chr'], row['Start'], interval)
        temp_h3k9me3 = information_extration(h3k9me3_df, current_sample, row['Chr'], row['Start'], interval)
        temp_rna = information_extration(rna_df, current_sample, row['Chr'], row['Start'], interval)
        
        
        if sum(temp_ctcf) >= 0.1 and sum(temp_h3k9me3) >=0.1 and sum(temp_rna) >= 0.1 and\
            len(temp_ctcf) == 1001 and len(temp_h3k9me3) == 1001 and len(temp_rna) == 1001:
        
        
            list_dnas += [[temp_dnas]]
            list_ctcf += [temp_ctcf]
            list_h3k9me3 += [temp_h3k9me3]
            list_rna += [temp_rna]
            
    
            list_name += [f'{row["Chr"]}_{row["Start"]}_{current_sample}']
            count += 1
            if count%500== 0:
                print(f'Train Finished {count}')
    
    
        if count >= ((samples_size/0.9) * 0.1):
            break
    
    
    
    result = pd.DataFrame(list_dnas, index = list_name).transpose()
    result.to_csv("val_DNase.csv", index=False)
    
    result = pd.DataFrame(list_ctcf, index = list_name).transpose()
    result.to_csv("val_ctcf.csv", index=False)
    
    
    result = pd.DataFrame(list_h3k9me3, index = list_name).transpose()
    result.to_csv("val_h3k9me3.csv", index=False)
    
    
    result = pd.DataFrame(list_rna, index = list_name).transpose()
    result.to_csv("val_rna.csv", index=False)
    
    
    
    
    #processing test samples
    random_df =df_test.sample(n = int((samples_size/0.9) * 2.0))
    list_ctcf = []
    list_h3k9me3 = []
    list_rna = []
    list_dnas = []
    list_name = []
    
    count = 0
    for index,row in random_df.iterrows():
        sample_index = random.randrange(0, num_sample - 1)   # skip first three line
        
        
        
        current_sample = samples[sample_index]
        
        temp_dnas = row[current_sample]
        temp_ctcf = information_extration(ctcf_df, current_sample, row['Chr'], row['Start'], interval)
        temp_h3k9me3 = information_extration(h3k9me3_df, current_sample, row['Chr'], row['Start'], interval)
        temp_rna = information_extration(rna_df, current_sample, row['Chr'], row['Start'], interval)
        
        
        if sum(temp_ctcf) >= 0.1 and sum(temp_h3k9me3) >=0.1 and sum(temp_rna) >= 0.1 and\
            len(temp_ctcf) == 1001 and len(temp_h3k9me3) == 1001 and len(temp_rna) == 1001:
        
            
        
        
            list_dnas += [[temp_dnas]]
            list_ctcf += [temp_ctcf]
            list_h3k9me3 += [temp_h3k9me3]
            list_rna += [temp_rna]
            
    
            list_name += [f'{row["Chr"]}_{row["Start"]}_{current_sample}']
            count += 1
            if count%500== 0:
                print(f'Train Finished {count}')
    
    
        if count >= ((samples_size/0.9) * 0.1):
            break
    
    
    
    
    result = pd.DataFrame(list_dnas, index = list_name).transpose()
    result.to_csv("test_DNase.csv", index=False)
    
    result = pd.DataFrame(list_ctcf, index = list_name).transpose()
    result.to_csv("test_ctcf.csv", index=False)
    
    
    result = pd.DataFrame(list_h3k9me3, index = list_name).transpose()
    result.to_csv("test_h3k9me3.csv", index=False)
    
    
    result = pd.DataFrame(list_rna, index = list_name).transpose()
    result.to_csv("test_rna.csv", index=False)
    
    
    
    
    
    
    
    
def information_extration(df, sample_name, chr_name, start_position, interval):
    
    df = df[df['Chr'] == chr_name]
    df = df[(df['Start'] >= (start_position - interval)) & (df['Start'] <= (start_position + interval))]

    
    result = df[sample_name].tolist()
     

    return result
    
    
    
    
    
    
    
    
"""
ctcf_df, h3k9me3_df, dnase_df, rna_df = read_df("ctcf.bin","h3k9me3.bin","dnase.bin","rna.bin")
    
    
sample_split(ctcf_df, h3k9me3_df, dnase_df, rna_df)


    


    
ctcf_df = pd.read_csv('outsample_ctcf.bin')

h3k9me3_df = pd.read_csv("outsample_h3k9me3.bin")

dnase_df = pd.read_csv("outsample_DNase.bin")

rna_df = pd.read_csv("outsample_rna.bin")
    
    
out_sample_process(ctcf_df, h3k9me3_df, dnase_df, rna_df, 10000, interval = 500000, random_seed = 123) 
    
    
    
ctcf_df = pd.read_csv('insample_ctcf.bin')

h3k9me3_df = pd.read_csv("insample_h3k9me3.bin")

dnase_df = pd.read_csv("insample_DNase.bin")

rna_df = pd.read_csv("insample_rna.bin")


random_sample_generation(ctcf_df, h3k9me3_df, dnase_df, rna_df, samples_size, interval = 500000, random_seed = 123)




"""


















    
    
    
    
    
    
    
    
    
    
    