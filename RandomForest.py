import pandas as pd
import numpy as np
import sys
from sklearn import datasets

from scipy import stats as stats

#for i in range (1,len(sys.argv),2):
#      if sys.argv[i] == "-df":
#        DF = sys.argv[i+1]
#df = pd.read_csv(DF, sep='\t',header=0) 


iris = datasets.load_iris()
x = iris.data[:,[2,3]]
y = iris.target

from sklearn.cross_validation import train_test_split
X_train, X_test, Y_train, Y_test = train_test_split(x,y,test_size=0.3, random_state=0)

from sklearn.ensemble import RandomForestClassifier 
forest = RandomForestClassifier(criterion='entropy',n_estimators=10, random_state=1, n_jobs=2)
forest = forest.fit(X_train, Y_train)
plot_decision_regions(X_combined, Y_combined, classifier=forest, test_idx=range(105,150))


#forest = forest.fit(train_data[],train_data[])