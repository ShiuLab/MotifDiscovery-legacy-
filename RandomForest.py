import pandas as pd
import numpy as np
import sys
from sklearn.ensemble import RandomForestClassifier 
from scipy import stats as stats

for i in range (1,len(sys.argv),2):
      if sys.argv[i] == "-df":
        DF = sys.argv[i+1]

df = pd.read_csv(DF, sep='\t',header=0) 

forest = RandomForestClassifier(n_estimators=100)
#forest = forest.fit(train_data[],train_data[])