import pandas as pd
import numpy as np
import os
    
def get_nid():
    nid={}
    index=0
    for i in "ATGC":
        for j in "ATGC":
            nid[i+j]=index
            index+=1 
    return nid
        
def chr_col_preprocessing_data(data):
    chr_data=np.array(data["chr"])
    out=[[],[],[]]
    for i in range(len(chr_data)):
        x=chr_data[i][1:]
        s=x.split(":")
        out[0].append(s[0][3:])
        s=s[1].split("-")
        # print(s)
        out[1].append(int(s[0]))
        out[2].append(int(s[1])+100)
    data["chr"]=out[0]
    data["start"]=out[1]
    data["end"]=out[2]
    
    return data

def seq_to_one_hot_data(data):
    nid=get_nid()
    seq_data=data["sequence"]
    one_hot=[np.zeros(len(seq_data)) for z in range(16)]
    for i in range(len(seq_data)):
        seq=seq_data.iloc[i]
        if seq in nid.keys():
            one_hot[nid[seq]][i]+=1
    index=0
    for i in nid.keys():
        data[i]=one_hot[index]
        index+=1
    
    return data.drop(["sequence"],axis=1,inplace=True)

def chr_col_preprocessing(data):
    chr_data=np.array(data["chr"])
    out=[[],[],[]]
    for i in range(len(chr_data)):
        x=chr_data[i][1:]
        s=x.split(":")
        out[0].append(s[0][3:])
        s=s[1].split("-")
        out[1].append(int(s[0]))
        out[2].append(int(s[1]))
    data["chr"]=out[0]
    data["start"]=out[1]
    data["end"]=out[2]
    
    return data

def seq_to_one_hot(data):
    nid=get_nid()
    seq_data=data["sequence"]
    one_hot=[np.zeros(len(seq_data)) for z in range(16)]
    for i in range(len(seq_data)):
        seq=seq_data.iloc[i]
        for j in range(len(seq)-1):
            if seq[j:j+2] in nid.keys():
                one_hot[nid[seq[j:j+2]]][i]+=1
    index=0
    for i in nid.keys():
        data[i]=one_hot[index]
        index+=1
    
    return data.drop(["sequence"],axis=1,inplace=True)


def generate_di_nucleotide_frequency_file(file_dir):

    os.chdir(file_dir)
    
    data_temp=pd.read_csv("output.fasta",header=None)
    data=pd.DataFrame()
    data["chr"]=data_temp.iloc[0::2]
    data["sequence"]=np.array(data_temp.iloc[1::2][0])
    
    chr_col_preprocessing_data(data)
    seq_to_one_hot_data(data)
    
    data_cpg_temp=pd.read_csv("/project/ReadStatistics-data/output_cpg.fasta",header=None)
    data_cpg=pd.DataFrame()
    data_cpg["chr"]=data_cpg_temp.iloc[0::2]
    data_cpg["sequence"]=np.array(data_cpg_temp.iloc[1::2][0])
    
    chr_col_preprocessing(data_cpg)
    seq_to_one_hot(data_cpg)
    data_cpg.reset_index(inplace=True)
    
    
    overlap_cpg_dinuc=pd.read_csv("output_reads.bed",delimiter="\t",header=None)
    overlap_cpg_dinuc[1]=overlap_cpg_dinuc[1]
    overlap_cpg_dinuc[2]=overlap_cpg_dinuc[2]+100
    overlap_cpg_dinuc[0]=overlap_cpg_dinuc[0].str.lstrip("chr")
    overlap_cpg_dinuc[3]=overlap_cpg_dinuc[3].str.lstrip("chr")
    overlap_cpg_dinuc.rename({0: 'chr', 1: 'start',2:"end"}, axis=1, inplace=True)
    overlap_cpg_dinuc = pd.merge(overlap_cpg_dinuc,data, on=['chr','start','end'])
    
    overlap_cpg_dinuc.drop(["start","end"],axis=1,inplace=True)
    df_new = overlap_cpg_dinuc.groupby([3,4,5],as_index=False).sum(numeric_only=True)
    df_new.rename({3: 'chr', 4: 'start', 5:"end"}, axis=1, inplace=True)
    
    df_final=pd.merge(data_cpg, df_new, on=['chr','start','end'])
    df_final.drop(["index"],axis=1,inplace=True)
    
    df_final.to_csv("FinalFile.csv",index=False)
    
    print("Done! :",i)
    