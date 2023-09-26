from __init__ import *

file_names = sys.argv[2:]
combined_methylation_outputs=sys.argv[1]

summary_file =  pd.DataFrame(columns=['chr','start','end'])
for bed_file in file_names:
    #print(bed_file)

    data = pd.read_csv(bed_file, sep="\t",low_memory=False)
    data.rename({"methylation_level_threshold_adjusted":bed_file.split("/")[-2].split(".")[0]+"_methylation_level_threshold_adjusted"   },axis=1,inplace=True)
    data.rename({"methylation_level":bed_file.split("/")[-2].split(".")[0]+"_methylation_level"  },axis=1,inplace=True)
    
    data['chr']=data['chr'].astype(np.str_)
    data['start']=data['start'].astype(np.int8)
    data['end']=data['end'].astype(np.int8)

    summary_file = pd.merge(data,summary_file,on=['chr','start','end'], how="outer")

summary_file.to_csv(combined_methylation_outputs,index=False,sep="\t")
