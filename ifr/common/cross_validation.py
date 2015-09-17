'''
Created on Sep 6, 2013
@author: Stephen O'Hara

Module with functions helpful to generate cross-validation
splits of data.

This module is a refactoring and extension of code previously
found in the svo_util.data_mgmt module.
'''
import scipy as sp

#======================================
# Top-level cross-validation generators
#======================================

def CV_Generator(n, func, *args, **kwargs):
    '''
    The function creates an iterator (generator) that yields
    (train_idxs, test_idxs) for use in cross-validating machine
    learning methods. The generator produced is intended to be
    compatible for use with sklearn functions that take CV objects.
    
    For example::
    
      S = sp.array([1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8])
      f = subject_partition_idxs
      cvg = CV_Generator(5, f, S, frac=0.25)
      for (train_idxs,test_idxs) in cvg:
          print "Train idxs: %s"%str(train_idxs)
          print "Test idxs: %s"%str(test_idxs)
    
    @param n: The number of splits to generate
    @param func: The splitting function that returns (train_idxs, test_idxs, _aux_info),
    such as subject_partition_idxs
    @param args: The arguments required by func
    @param kwargs: The keyword arguments required by func
    @return: A cross-validation index generator.
    '''
    count = 0
    while count < n:
        (t1, t2, _) = func(*args,**kwargs)
        yield (t1, t2)
        count += 1
        

def StratifiedShuffleSplitSubjects(n, L, S, frac):
    '''
    This function produces a generator that yields
    train/test index sets based on the subject_partition_proportional()
    splitting methodology. This pattern should be compatible with scikits.learn (sklearn)
    functions that take a cross validation iterator as input.
    @note: This is just a convenience wrapper to using CV_Generator(...)
    @param n: The number of partitions to generate
    @param L: The vector of labels, one integer per sample
    @param S: The vector of subject ids, one integer per sample
    @param frac: The fraction of the data used for testing.
    @return: A generator that can be used to iterate over the splits
    '''
    f = subject_partition_proportional_idxs
    cvg = CV_Generator(n, f, L, S, frac=frac)
    return cvg

#======================================
# Top-level functions for applying the
# selected train/test indexes to split
# the Data (D) and Labels (L) matrices
#======================================

def subject_partition(D, L, S, frac=None, subject_id=None):
    '''
    Partitions data to leave one or more subjects out. This is useful when
    you have data with subject overlap, and you wish to cross validate
    by leaving out all samples for a given subject id.
    @param D: The data matrix, samples are in rows. Features in columns.
    @param L: A vector of target labels.
    @param S: A vector of subject ids, which indicate which subject the
    data was sampled from.
    @param frac: Optional. Specify a fraction, such as 0.25, which indicates
    the fraction of subjects to randomly leave out. If none, then only one
    subject id will be selected.
    @param subject_id: Optional. Set this parameter to an integer to
    specify the subject_id to use in the loso partition.
    If None, then a random subject will be chosen. NOTE: Only used if frac is None. 
    @return: ( train_data, train_labels, test_data, test_labels, test_ids)
    '''
    rc = subject_partition_idxs(S, frac=frac, subject_id=subject_id)
    return _partition_data(D, L, rc)

def subject_partition_proportional(D, L, S, frac=None):
    '''
    This is similar to subject-based splitting, but it first
    stratifies the data into the class labels (unique entries in L),
    and performs subject partitioning in each stratification.
    This allows the user to have train/test splits that have no
    subject id overlap while also having approximately the same
    proportion of class labels.
    @return: ( train_data, train_labels, test_data, test_labels, test_ids)
    '''
    rc = subject_partition_proportional_idxs(L, S, frac=frac)
    return _partition_data(D, L, rc)

def random_partition(D, L, pctTrain=0.75, proportional_labels=True):
    '''
    Partitions Data matrix and associated Label vector into
    randomized Train and Test partitions.
    @param D: The data matrix, samples in rows, as a numpy array.
    @param L: The label vector corresponding to D, as a numpy array.
    @param pctTrain: The percentage of samples to include in training partition
    @param proportional_labels: If true, the pctTrain will be applied to each
    label. If false, pctTrain will be applied to the entire data set, ignoring labels.
    Only with large data sets having approximately equal labels should you not
    do proportional partitioning.
    @return: (Dtrain, Ltrain, Dtest, Ltest, p), where p is the permutation
    '''
    
    #note: for backwards compatibility with old code, pctTrain is used in
    # this top-level function whereas the others use frac (test fraction)
    # below we convert pctTrain to frac before calling supporting function.
    rc = random_partition_idxs(L, frac=(1-pctTrain), proportional_labels=proportional_labels)
    return _partition_data(D, L, rc)


#======================================
# Core functions that produce a single
# (train/test) split, where only the
# indexes are returned. It is left to
# the user (or a higher level function)
# to apply the indices to the data.
#======================================

def subject_partition_idxs(S, frac=None, subject_id=None):
    '''
    Creates a cross-validation partitioning by splitting the
    indexes sensitive to the subject labels (to prevent subject overlap).
    This generates 'leave one subject out' or 'leave p subjects out'
    @param S: The subject ids per sample. len(S) is the number of
    samples in the data set
    @param frac: If None, then leave-one-subject out is performed. Else,
    leave p subjects out where p is the fraction of the unique subject labels.
    @param subject_id: Ignored if frac is not None. This is the specific
    subject_id to be used for testing. If None, then a single random id is
    chosen uniformly over the set of unique subject ids.
    @return (train indexes, test indexes, test ids )
    '''
    unique_s = sp.unique(S)
    if not subject_id is None:
        assert subject_id in S, "Subject id: %d is not an element of S."%subject_id
        
    if frac is None:
        #leave one subject out
        if subject_id is None:
            p = sp.random.random_integers(0,len(unique_s)-1)  #pick a random index from [low,high]
            loso_id = unique_s[p]
        else:
            loso_id = subject_id
        test_ids = [loso_id]     
    else:
        #leave more than one subject out
        p = sp.random.permutation(range(len(unique_s)))
        Ntest = int(frac*len(unique_s))
        test_ids = sorted(unique_s[p[0:Ntest]])
        
    idxs = range(len(S))
    test_idxs = sorted( [ i for i in idxs if S[i] in test_ids])
    train_idxs = sorted( set(idxs) - set(test_idxs) )
    return (train_idxs, test_idxs, test_ids)

def subject_partition_proportional_idxs(L,S, frac=None):
    '''
    This is similar to subject-based splitting, but it first
    stratifies the data into the class labels (unique entries in L),
    and performs subject partitioning in each stratification.
    This allows the user to have train/test splits that have no
    subject id overlap while also having approximately the same
    proportion of class labels.
    @return: (train_idxs, test_idxs, test_ids)
    '''
    test_ids = []
    for label in sp.unique(L):
        Sx = S[L==label]  #subjects for a specific class, choose frac of these
        tmp = subject_partition_idxs(Sx, frac=frac)
        #returns indexes into Sx, so we'll use the subject ids instead
        ids = _get_test_subjects(Sx, tmp )
        test_ids += list(ids)
    idxs = range(len(S))
    test_idxs = sorted( [i for i in idxs if S[i] in test_ids] )
    train_idxs = sorted( set(idxs) - set(test_idxs))
    return (train_idxs, test_idxs, test_ids)


def random_partition_idxs(L, frac=0.25, proportional_labels=True):
    '''
    Generates indexes to partition data for train/test splits. This
    function can generate random partitioning (ShuffleSplit) or
    random partitioning where the labels are kept proportional (Stratified Shuffle).
    @param L: The label vector as a numpy array or list. This is used to get the
    number of samples in the data len(L), and also for data stratification along
    class labels if proportional_labels is True.
    @param frac: The percentage of samples to include in test partition
    @param proportional_labels: If true, the frac will be applied to each
    label. If false, frac will be applied to the entire data set, ignoring labels.
    Only with large data sets having approximately equal labels should you NOT
    do proportional partitioning.
    @return: (train_idxs, test_idxs, p), where p is the random permutation or
    list of permutations (one per class) if proportional is True.
    '''
    if not proportional_labels:
        N = len(L)
        p = sp.random.permutation(range(N))
        Ntest = int(frac*N)
        test_idxs  = sorted( p[0:Ntest] )
        train_idxs = sorted( p[Ntest:]  )
    else:
        test_idxs = []
        p = []
        for label in sp.unique(L):  #this also sorts the unique labels
            Lx = sp.flatnonzero(L==label) #indexes in L where L == label
            N = len(Lx)
            px = sp.random.permutation(range(N))
            Ntest = int(frac*N)
            test_idxs += list( Lx[px[0:Ntest]])
            p.append(px)
        test_idxs = sorted(test_idxs)
        train_idxs = sorted( set(range(len(L))) - set(test_idxs))
            
    return (train_idxs, test_idxs, p)


#======================================
#Private helper functions to reduce
# code duplication in above code.
#======================================
def _partition_data(D,L,partition_idxs):
    (train_idxs, test_idxs, aux_info) = partition_idxs
    train_data = D[train_idxs,:]
    train_labels = L[train_idxs]
    test_data = D[ test_idxs, : ]
    test_labels = L[ test_idxs ]
    return (train_data, train_labels, test_data, test_labels, aux_info)

def _get_test_subjects(S, partition_idxs):
    test_idxs = partition_idxs[1]
    return sp.unique( sp.array(S)[test_idxs] )
