from __init__ import *
from MachineLearning import *

model_dir=sys.argv[1]
finalfile_dir=sys.argv[2]
methylation_annotation=sys.argv[3]


ML=MachineLearning()
ML.train_save(model_dir,finalfile_dir,methylation_annotation)
