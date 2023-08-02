from __init__ import *

class MachineLearning:
    
    def __init__(self):
        
        self.originalclass = []
        self.predictedclass = []
        
    def classification_report_with_accuracy_score(self,y_true, y_pred):
        
        self.originalclass.extend(y_true)
        self.predictedclass.extend(y_pred)
        
        return accuracy_score(y_true, y_pred)

    def normalize(self,x):
        new_x=[]

        for i in range(len(x)):
            x_new=[]
            for j in range(16,32):
                a=np.sum(x[i][1:17])
                b=np.sum(x[i][17:33])
                if x[i][j-16]==0 or b==0:
                    x_new.append(0.0)
                else:
                    x_new.append(x[i][j]*a/( x[i][j-16]*b))
            new_x.append(x_new)
        return new_x

    def confidence_score(self,y_true,y_pred_prob):

        keys=[0,1]
        dict={}
        i=0
        for key in keys:
            dict[key]=i
            i+=1

        cs=0.0
        for i in range(len(y_true)):
            cs+=y_pred_prob[i][dict[y_true[i]]]


        return cs/len(y_true)
        
    def data_extraction_for_ML(self,final_file):
        
        x=final_file
        y=np.array(x[cancer_type])
        y=np.array(y).reshape(len(y)  ,)
        y_1=np.where(y<0.5)[0]
        y_2=np.where(y>0.5)[0]
        index=np.random.choice(y_1,len(y_2), replace=False)
        index=np.concatenate((index, y_2), axis=0)
        x_32=x[[ 'AA_x', 'AT_x', 'AG_x', 'AC_x',
           'TA_x', 'TT_x', 'TG_x', 'TC_x', 'GA_x', 'GT_x', 'GG_x', 'GC_x', 'CA_x',
           'CT_x', 'CG_x', 'CC_x', 'AA_y', 'AT_y', 'AG_y', 'AC_y', 'TA_y', 'TT_y',
           'TG_y', 'TC_y', 'GA_y', 'GT_y', 'GG_y', 'GC_y', 'CA_y', 'CT_y', 'CG_y',
           'CC_y']]
        y_32=y
        x_32_norm=np.array(self.normalize(np.array(x_32)))
        x, y = shuffle(x_32_norm[index], y[index], random_state=0)
        y[y<0.5]=0.0
        y[y>0.5]=1.0

        return x,y
        
    def random_forest_model(self,x,y):
        
        x=StandardScaler().fit_transform(x)
        clf_rf = RandomForestClassifier(class_weight='balanced', max_depth=20 ,n_estimators=800, criterion='entropy')
   
        outer_cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=1)
        nested_score = cross_val_score(clf_rf, X=x, y=y, cv=outer_cv, scoring=make_scorer(self.classification_report_with_accuracy_score))

        print(classification_report(self.originalclass, self.predictedclass))
    
        return clf_rf 
    
    def data_extraction_for_ML(self,x):
        x_32=x[[ 'AA_x', 'AT_x', 'AG_x', 'AC_x',
           'TA_x', 'TT_x', 'TG_x', 'TC_x', 'GA_x', 'GT_x', 'GG_x', 'GC_x', 'CA_x',
           'CT_x', 'CG_x', 'CC_x', 'AA_y', 'AT_y', 'AG_y', 'AC_y', 'TA_y', 'TT_y',
           'TG_y', 'TC_y', 'GA_y', 'GT_y', 'GG_y', 'GC_y', 'CA_y', 'CT_y', 'CG_y',
           'CC_y']]
        x_32_norm=np.array(self.normalize(np.array(x_32)))

        x_norm =StandardScaler().fit_transform(x_32_norm)
   
        return x_norm
    def predict_save(self,model_dir,finalfile_dir,methylation_pred_dir):
        model = load(model_dir) 
        data=pd.read_csv(finalfile_dir,low_memory=False)
        x=self.data_extraction_for_ML(data)
        data['methylation_level']=model.predict_proba(x)[:,1]
        data[['chr','start','end','methylation_level']].to_csv(methylation_pred_dir,index=False)
