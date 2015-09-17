'''
Created on Dec 17, 2012
@author: Stephen O'Hara

This module defines a standard interface and implementation
for a variety of machine learning (ml) algorithms for use
as "ml_engines" for other modules, such as feature learning.
'''
import ifr

import scipy as sp
try:
    from sklearn import svm as svm
    from sklearn import linear_model as lm
except:
    print "Error: sklearn package not available."
    print "Machine learning functions will not work."

def svm_engine(Train, Test, perm=None, common_scale=True, verbose=False, no_normalization=False,
               loss='l2', penalty='l1', dual=False, C=0.5, tol=0.0001, **kwargs):
    '''
    svm randomized feature extraction engine
    Internal function that performs the "engine" of the sparse linear svc
    fitting of the data.
    @param Train: A tuple (D,L) providing the data and labels for the training
    data, for example, as loaded by the idr.load_flu_mat function.
    @param Test: A tuple (D,L) for the test data.
    @param perm: The random factor/column permutation of the data. None means
    no permutation. Otherwise this is a permuted list of column indexes. You
    can also use this parameter to perform non-permuted subset selection. Just
    provide a list of column indexes to include.
    @param common_scale: If True, then the test data will be normalized using the centers/scale
    of the training data (i.e. "common scale"). This is a common thing to do in machine learning --
    project the test samples onto the normalized coordinates of the training data before applying
    the model. This breaks, however, when training on H3N2 microarray data and testing on H1N1
    because the collection procedures result in very different scale expression data. It is better
    in this case to normalize each separately, and use the self-normalized test vectors with
    the trained model.
    @param verbose: If True, incremental output to the console will be generated. Otherwise,
    the results will be returned silently.
    @param no_normalization: If True, then the built-in normalizing done by this function will
    be skipped. This is appropriate if the data being provided is already normalized.
    @return: A tuple (test_accuracy, factors, clf, train_accuracy) where test_accuracy is
    the accuracy on the test data, factors are those columns with non-zero coefficients used
    in the classifier, clf is the trained classifier, and train_accuracy is the accuracy
    on the training data.
    
    @note: Other parameters, loss, penalty, dual, C, tol, are input parameters to the scikit-learn
    SVC classifier. Please see the documentation on the LinearSVC for definitions.
    '''
    (D,L) = Train
    (D2, L2) = Test
    
    if not perm is None:
        D = D[:,perm]
        D2 = D2[:,perm]
    
    if no_normalization:
        Xtrn = D
        Xtst = D2
    else:
        #normalize the training data, mean-center unit-deviation
        (Xtrn,means,stdvs) = ifr.normalize_data(D)
        
        if common_scale:
            #use the training normalization means/stdvs to apply to testing data
            (Xtst,_,_) = ifr.normalize_data(D2, means, stdvs)
        else:
            #use the test statistics to normalize the test data separately from the training statistics
            (Xtst,_,_) = ifr.normalize_data(D2)
       
    
    #parameter C=1.0 determined by a parameter search on the training data
    # using leave-one-subject-out...see svm_find_best_param function
    # seems to work well even for very different sizes of the number of columns in the data
    clf = svm.LinearSVC(loss=loss, penalty=penalty, dual=dual, C=C, tol=tol, **kwargs) #, class_weight='auto')   
    clf.fit(Xtrn,L)
    
    #print "Non-zero coefficients have values:"
    #print clf.raw_coef_[ sp.nonzero(clf.raw_coef_)]
    
    if verbose: print "Predicting training values to measure training error..."
    s1 = clf.score(Xtrn,L)
    
    if verbose:
        print "%d correct predictions on training samples"%( s1*len(L) )
        print "%.3f fraction correct of %d samples"%( s1, len(L) )
    
    if verbose: print "Predicting test values"
    s2 = clf.score(Xtst, L2)
    preds = clf.predict(Xtst)
    if verbose:
        print "%d correct predictions on testing samples"%(s2*len(L2))
        print "%.3f fraction correct of %d samples"%( s2, len(L2))
    
    x = list(clf.coef_[0,:])
    factors = sp.nonzero(x)[0]
    if not perm is None:
        if type(perm) is list:
            perm = sp.array(perm) #make a scipy array
        factors = perm[factors]
    if verbose: print "%d gene indexes used in model:"%len(factors)
    if verbose: print factors
    
    return (s2, factors, clf, s1, preds)


def lr_engine(Train, Test, perm=None, common_scale=True, verbose=False, no_normalization=False,
              penalty='l1', intercept_scaling=10, C=1.0, tol=0.00001, **kwargs):
    '''
    logistic regression (lr) classification engine
    Internal function that performs the "engine" of the sparse linear regression
    fitting of the data.
    @param Train: A tuple (D,L) providing the data and labels for the training
    data, for example, as loaded by the idr.load_flu_mat function.
    @param Test: A tuple (D,L) for the test data.
    @param perm: The random factor/column permutation of the data. None means
    no permutation. Otherwise this is a permuted list of column indexes. You
    can also use this parameter to perform non-permuted subset selection. Just
    provide a list of column indexes to include.
    @param common_scale: If True, then the test data will be normalized using the centers/scale
    of the training data (i.e. "common scale"). This is a common thing to do in machine learning --
    project the test samples onto the normalized coordinates of the training data before applying
    the model. This breaks, however, when training on H3N2 microarray data and testing on H1N1
    because the collection procedures result in very different scale expression data. It is better
    in this case to normalize each separately, and use the self-normalized test vectors with
    the trained model.
    @param verbose: If True, incremental output to the console will be generated. Otherwise,
    the results will be returned silently.
    @return: A tuple (test_accuracy, factors, clf, train_accuracy) where test_accuracy is
    the accuracy on the test data, factors are those columns with non-zero coefficients used
    in the classifier, clf is the trained classifier, and train_accuracy is the accuracy
    on the training data.
    
    @note: Other parameters are input parameters to the scikit-learn logistic regression
    classifier (with L1 regularization). Please see the documentation on the LinearSVC for definitions.
    '''
    (D,L) = Train
    (D2, L2) = Test
    
    if not perm is None:
        D = D[:,perm]
        D2 = D2[:,perm]
    
    if no_normalization:
        Xtrn = D
        Xtst = D2
    else:
        #normalize the training data, mean-center unit-deviation
        (Xtrn,means,stdvs) = ifr.normalize_data(D)
        
        if common_scale:
            #use the training normalization means/stdvs to apply to testing data
            (Xtst,_,_) = ifr.normalize_data(D2, means, stdvs)
        else:
            #use the test statistics to normalize the test data separately from the training statistics
            (Xtst,_,_) = ifr.normalize_data(D2)
         
    clf = lm.LogisticRegression(penalty=penalty, C=C, fit_intercept=True, 
                                intercept_scaling=intercept_scaling, tol=tol, **kwargs)
    
    clf.fit(Xtrn,L)
    
    #print "Non-zero coefficients have values:"
    #print clf.raw_coef_[ sp.nonzero(clf.raw_coef_)]
    
    if verbose: print "Predicting training values to measure training error..."
    s1 = clf.score(Xtrn,L)
    
    if verbose:
        print "%d correct predictions on training samples"%( s1*len(L) )
        print "%.3f fraction correct of %d samples"%( s1, len(L) )
    
    if verbose: print "Predicting test values"
    s2 = clf.score(Xtst, L2)
    preds = clf.predict(Xtst)
    if verbose:
        print "%d correct predictions on testing samples"%(s2*len(L2))
        print "%.3f fraction correct of %d samples"%( s2, len(L2))
    
    x = list(clf.coef_[0,:])
    factors = sp.nonzero(x)[0]
    if not perm is None:
        if type(perm) is list:
            perm = sp.array(perm) #make a scipy array
        factors = perm[factors]
    if verbose: print "%d gene indexes used in model:"%len(factors)
    if verbose: print factors
    
    return (s2, factors, clf, s1, preds)

ML_ENGINE_SVM = svm_engine
ML_ENGINE_LR  = lr_engine


if __name__ == '__main__':
    pass