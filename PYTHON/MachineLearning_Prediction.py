from __init__ import *
from MachineLearning import *
import sys

model_dir=sys.argv[1]
finalfile_dir=sys.argv[2]
methylation_pred_dir=sys.argv[3]

ML=MachineLearning()
ML.predict_save(model_dir,finalfile_dir,methylation_pred_dir)