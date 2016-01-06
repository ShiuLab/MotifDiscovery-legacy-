"""
Testing TPOT [Tree-based Pipeline Optimization Tool] built by Randy Olson 
(http://www.randalolson.com/2015/11/15/introducing-tpot-the-data-science-assistant/)
"""

from tpot import TPOT
from sklearn.datasets import load_digits  
from sklearn.cross_validation import train_test_split  
  
digits = load_digits()  
X_train, X_test, y_train, y_test = train_test_split(digits.data, digits.target,  
                                                    train_size=0.75)  
  
tpot = TPOT(generations=5, verbosity=2)  
tpot.fit(X_train, y_train)  
print(tpot.score(X_train, y_train, X_test, y_test))
tpot.export('tpot_exported_pipeline.py')