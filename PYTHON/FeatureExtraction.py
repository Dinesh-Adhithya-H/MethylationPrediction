from __init__ import *
import sys

class FeatureExtraction:
    
    def __init__(self,dinucletide_data_dir,cpg_island_data_dir,overlap_cpg_dinucleotide_dir):
        
        self.dinucletide_data=pd.read_csv(dinucletide_data_dir,header=None)
        self.cpg_island_data=pd.read_csv(cpg_island_data_dir,header=None)
        self.overlap_cpg_dinucleotide=pd.read_csv(overlap_cpg_dinucleotide_dir,delimiter="\t",header=None)
        
    def process_files(self):
        data=pd.DataFrame()
        data["chr"]=self.dinucletide_data.iloc[0::2]
        data["sequence"]=np.array(self.dinucletide_data.iloc[1::2][0])
        
        self.chr_col_preprocessing_data(data)
        self.seq_to_one_hot_data(data)
        
        print("DINUCLEOTIDE DATA PROCESS DONE!")
        
        data_cpg=pd.DataFrame()
        data_cpg["chr"]=self.cpg_island_data.iloc[0::2]
        data_cpg["sequence"]=np.array(self.cpg_island_data.iloc[1::2][0])
        
        self.chr_col_preprocessing(data_cpg)
        self.seq_to_one_hot(data_cpg)
        
        print("CpG ISLAND DATA PROCESS DONE!")
        
        df_final = self.combine_data_and_cpg_dinuc(self.overlap_cpg_dinucleotide,data,data_cpg)
        
        print("DINUCLEOTIDE DATA AND CpG ISLAND DATA COMBINATION PROCESS DONE!")
        

        return df_final

    def chr_col_preprocessing_data(self,data):
        chr_data=np.array(data["chr"])
        out=[[],[],[]]
        for i in range(len(chr_data)):
            x=chr_data[i][1:]
            s=x.split(":")
            out[0].append(s[0][3:])
            s=s[1].split("-")
            out[1].append(int(s[0]))
            out[2].append(int(s[1])+100)
        data["chr"]=out[0]
        data["start"]=out[1]
        data["end"]=out[2]

        return data

    def seq_to_one_hot_data(self,data):

        nid={}
        index=0
        for i in "ATGC":
            for j in "ATGC":
                nid[i+j]=index
                index+=1     


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

    def chr_col_preprocessing(self,data):
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

    def seq_to_one_hot(self,data):

        nid={}
        index=0
        for i in "ATGC":
            for j in "ATGC":
                nid[i+j]=index
                index+=1     

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


    def rename(self,l):
        l_new=[]
        for i in l:
            if i=="X":
                l_new.append(23)
            elif i=="Y":
                l_new.append(24)
            else:
                l_new.append(int(i))
        return l_new

    def combine_data_and_cpg_dinuc(self,overlap_cpg_dinuc,data,data_cpg):
        
        
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
        
        return df_final



    def methylation_combine_cpg(self,df_final,methyl):
    
        df_final.chr=self.rename(df_final.chr)
        methyl.chr=self.rename(methyl.chr)
        methyl=methyl[['chr','start','end', self.cancer_name]]
        final_file=pd.merge(methyl,df_final,on=['chr','start','end'])

        return final_file
    
dinucletide_data_dir = sys.argv[1] # dinucletide_data_dir
cpg_island_data_dir = sys.argv[2] # cpg_island_data_dir
overlap_cpg_dinucleotide_dir = sys.argv[3] # overlap_cpg_dinucleotide_dir
output_dir_name=sys.argv[4]

Feature_Extraction=FeatureExtraction(dinucletide_data_dir,cpg_island_data_dir,overlap_cpg_dinucleotide_dir)
final_file=Feature_Extraction.process_files()
final_file.to_csv(output_dir_name+"final_file.csv")