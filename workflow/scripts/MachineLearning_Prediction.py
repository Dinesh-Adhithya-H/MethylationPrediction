from __init__ import *
from MachineLearning import *

model_dir=sys.argv[1]
finalfile_dir=sys.argv[2]
methylation_pred_dir=sys.argv[3]
ratio=sys.argv[4]

ML=MachineLearning()
ML.predict_save(model_dir,finalfile_dir,methylation_pred_dir,ratio)
