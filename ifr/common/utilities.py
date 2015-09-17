'''
Created on Sep 25, 2012
@author: Stephen O'Hara

Misc. utility functions that
are used throughout the code base.
'''
import sys
import scipy as sp
import scipy.stats as stats
import urllib
from StringIO import StringIO
import tokenize

def fetch_data(url, params, method='POST'):
    '''
    Convenient wrapper around urllib functions to
    post/get a url to a web server and get the results.
    '''
    params = urllib.urlencode(params)
    if method=='POST':
        f = urllib.urlopen(url, params)
    else:
        f = urllib.urlopen(url+"?"+params)
    return (f.readlines(), f.code)

def _reporthook(a,b,c): 
    # ',' at the end of the line is important!
    print "% 3.1f%% of %d bytes\r" % (min(100, float(a * b) / c * 100), c),
    sys.stdout.flush()
    
def download_file(url, dest_file):
    '''
    Downloads a file given by a url (such as one of the zipfiles
    for the kent ridge biomedical repository), and stores the
    file to the destination directory.
    @param url: The full url to the file to be downloaded. ifr.LUNG_CANCER_URL,
    ifr.PROSTATE_URL, and ifr.BCELL_LYMPHOMA_URL, are constants that point
    to the appropriate file.
    @param dest_file: The destination filename, including full path.
    @return: (filename, headers)
    '''
    rc = urllib.urlretrieve(url, dest_file, reporthook=_reporthook)
    return rc

def get_official_name(gene, entrezgene=False):
    '''
    For a given gene symbol, like ARTS-1, return the official
    gene name, such as ERAP1. This is done by a web-query to ncbi,
    and parsing the output. Thus, it is relatively slow and should only
    be used when necessary.
    @param gene: A string like 'OAS1', or entrez-gene id like '4010'
    if entrezgene=True
    @param entrezgene: If True, then a different query format is used
    to resolve the gene symbol. If the input is an entrezgene id,
    then this flag should be True, else False.
    @return: A string, representing the official name of the gene.
    '''
    if entrezgene:
        url_pre = "http://www.ncbi.nlm.nih.gov/gene/?term="
        url_post = "&format=text"
    else:
        url_pre = "http://www.ncbi.nlm.nih.gov/gene/?term="
        url_post = "+AND+(%22Homo+sapiens%22%5BOrganism%5D)&format=text"
    #url_post = "%5Bsym%5D+AND+(%22Homo+sapiens%22%5BOrganism%5D)&format=text"
    qurl = "%s%s%s"%(url_pre,gene,url_post)
    rc = fetch_data(qurl, "")
    if rc[1] != 200:
        print "Error retrieving information for gene: %s"%gene
        return None
    
    lines = rc[0]  #text returned from http post
    tmp_line = None
    for x in range(len(lines)):
        if 'Official Symbol' in lines[x]:
            tmp_line = lines[x]
            break
        
    if not tmp_line is None:
        new_symb = tmp_line.split(": ")[1].split(" ")[0]
        return new_symb
    else:
        return None

def rolling_average(input_list, window_size=10):
    '''
    Computes the rolling averages of the input list
    over the window_size...essential the average of
    the list elements surrounding the current element
    '''
    output_list = []
    radius = window_size/2
    for idx, _num in enumerate(input_list):
        min_idx = max(0, idx-radius)
        max_idx = min(len(input_list), idx+radius+1)
        tmp = input_list[min_idx:max_idx]
        #print "DEBUG: idx=%d, span=(%d,%d)"%(idx, min_idx, max_idx)
        output_list.append( sp.mean(tmp))
    return output_list

def print_progress(cur, total):
    '''
    This function can be called in a processing loop
    to print out a progress indicator represented
    by up to 10 lines of dots, where each line 
    represents completion of 10% of the total iterations.
    @param cur: The current value (integer) of the iteration/count.
    @param total: The total number of iterations that must occur.
    '''
    one_line = 40 if total < 400 else round( total / 10.0 )
    one_dot = 1 if one_line / 40.0 < 1 else round( one_line / 40.0)    
    if (cur+1)%one_line == 0:
        print ' [%d]'%(cur+1)
    elif (cur+1)%one_dot == 0:
        print '.',
        sys.stdout.flush()    
    if cur+1 == total: print ""

def parse_number(s, fail=None):
    '''
    Tries to convert input object s into a float. If it
    can, the float representation is returned. If not,
    then either it returns None, or returns a user-specified
    value.
    @param s: Input object to be converted to float()
    @param fail: Returned when s can't be converted to float.
    @return: Either float(s) or the value specified by fail, 
    default None.    
    '''
    try:
        f = float(s)
    except ValueError:
        return fail
    return f


def smart_split(linein, sep=","):
    '''
    Works much like built-in split method of strings, but this version
    will honor quoted strings. Meaning that if the separator is inside a quoted string,
    then the separator is ignored in favor of keeping quoted strings as a single token.
    '''
    curtokens = []
    fieldvals = []
    prev_end = 0
    tkz = tokenize.generate_tokens(StringIO(linein).readline)
    for _, tokval, (_,colstart), (_,colend), _ in tkz:
        ws_delta = colstart - prev_end
        prev_end = colend
        
        if tokval.strip() == '': continue  #ignore whitespace
        
        if ''.join(tokval.split(sep)) == '':
            if len(curtokens) > 0:
                fieldvals.append( ("".join(curtokens)).strip() )
                curtokens = []
                continue
            else:
                continue
        
        if ws_delta > 0:
            tokval = (" "*ws_delta)+tokval
        
        if (tokval[0] == tokval[-1] == '"') or (tokval[0] == tokval[-1] == '\''):
            tokval = tokval[1:-1]
            
        curtokens.append(tokval)
    #at end of line, we will likely have curtokens and no separator encountered
    if len(curtokens) > 0: fieldvals.append( ("".join(curtokens)).strip() )
        
    return fieldvals


def constant_cols(D):
    '''
    Generate an array that lists those columns that have
    zero standard deviation in the input matrix.
    @param D: The input data matrix, features in columns
    @return: A list of those columns which have zero std dev
    '''
    s = D.std(axis=0)
    idxs = sorted( [ i for i in range(len(s)) if s[i]==0])
    return idxs

def indicator_vars(L):
    '''
    Transforms the list of unique class labels, L, into
    an indicator matrix. One row in the matrix for each label.
    One column for each unique class (alphabetically sorted),
    such that a 1 in the column indicates the class
    @return: tuple (IV, unique_labels) where IV is the indicator
    variable matrix and unique_labels is the ordered set of 
    labels that correspond with the indicator columns.
    '''
    unique_labels = sorted( list(set(L)) )
    P = len(unique_labels)
    N = len(L)
    IV = sp.zeros((N,P))
    for i in range(N):
        tmp_label = L[i]
        iv_column = unique_labels.index(tmp_label)
        IV[i,iv_column] = 1
        
    return IV.astype(int), unique_labels

def normalize_data(D, means=None, stdvs=None):
    '''
    Transforms data matrix D such that columns (features) are mean-centered
    and unit-deviation. Remember to transform your testing data using
    the means and stdvs computed from your training data!!!
    @param D: NxP data matrix, rows are samples, features are columns
    @param means: An array of P values, representing the mean value that
    should be subtracted from each feature. Default is None, indicating the
    means should be computed from the samples in D.
    @param stdvs: An array of P values, representing the std deviations of
    the P columns to be used to normalize the columns. Default is None, indicating
    the stdevs will be computed from the samples in D.
    @return (Dn, means, stdvs) where Dn is NxP, such that columns are mean-centered
    and have unit deviation. means is the list of P mean values, and stdvs is the 
    list of standard deviations for each column
    '''
    means = sp.mean(D,axis=0) if means is None else means
    stdvs = sp.std(D,axis=0) if stdvs is None else stdvs
    Dn = (D - means) / stdvs
    return (Dn, means, stdvs)


def replace_missing_values_with_col_means(D, sentinel=-100000):
    '''
    For a numeric data matrix, D, cells which have missing
    data will be replaces with the column average for the cell's column.
    @param D: A numpy/scipy array. Assumes rows are samples and columns are field values
    @param sentinel: This is the special number that denotes a missing value. This
    must have been added to the data matrix by another function.
    @return: Dhat, where missing values are now mean values for the column.
    '''
    D2 = D.copy()
    _rows, cols = D.shape
    for c in range(cols):
        S = D[:,c] #array slice
        S2 = D2[:,c]
        mask1 = (S == sentinel)
        if sum( mask1 ) > 0:
            mask2 = (S != sentinel)
            S2[mask1] = sp.mean( S[mask2] )

    return D2
    
def replace_missing_values_with_value(D, sentinel=-100000, val=0):
    '''
    For a numeric data matrix, D, cells which have missing
    data will be replaces with the constant value specified
    in the arguments.
    @param D: A numpy/scipy array. Assumes rows are samples and columns are field values
    @param sentinel: This is the special number that denotes a missing value. This
    must have been added to the data matrix by another function.
    @param val: This is the number to use in place of missing entries
    @return: Dhat, where missing values are now mean values for the column.
    '''
    D2 = D.copy()
    _rows, cols = D.shape
    for c in range(cols):
        S = D[:,c] #array slice
        S2 = D2[:,c]
        mask1 = (S == sentinel)
        S2[mask1] = val

    return D2    
    
def gappy_means_and_stdvs(X, sentinel=0):
    '''
    Computes the column means and standard deviations of a data matrix X, where some of the values
    are missing or unknown. This function simply ignores those entries marked
    with the sentinel value in the computation.
    @param X: Data matrix, the columns of which will be processed, rows are assumed to be multiple
    observations.
    @param sentinel: The numeric value that represents a missing entry in X.
    @return: (M, S, N), where the mean and standard deviation of a column are computed using
    only the non-missing entries. I.e., mean(Y) = (sum Yi)/Nhat where i are the non-missing indexes, and
    Nhat is the number of non-missing entries. M is the vector of means, S is the vector of std devs,
    and N has the number of non-missing values per column
    '''
    (_N,p) = X.shape
    M = sp.zeros(p)  #vector of mean values
    S = sp.zeros(p)  #vector of std devs
    N = sp.zeros(p)  #vector of N, the number of non-missing values per column
    for col in range(p):
        v = X[:,col]
        x = v[ v != sentinel ]
        M[col] = x.mean()
        S[col] = x.std()
        N[col] = len(x)
        
    return (M,S,N)

def ttest_features(D,L, return_mat =False):
    '''
    Performs univariate t-test on each feature (column) in D,
    to determine the significance in the difference between
    the values for two classes.
    @param D: The data. Rows are samples, Columns are features
    @param L: The label vector of integers. There should be only two unique
    values (two classes) in L.
    @param return_mat: If true, then there are two return values, one is the
    sorted list of features based on t-measure, the second is a matrix in the
    same index order as D, with the T-measure and p-value for each feature. If
    false (default), only the sorted list is returned.
    @return: (y,TT), where y is a rank-ordered list of tuples, with the first in the list
    being the feature with highest ttest score. [ (score,pval,feature_idx),
    (score2, pval2, feature_idx2), .... (scoreN, pvalN, feature_idxN)]; TT is the
    matrix, in the same index order as the samples in D and L, which has in the
    first column, the t-measure, and the p-value in the second.
    '''
    #divide D into two matrices, one per class
    classes = sp.unique(L)
    assert len(classes) == 2
    classA = classes[0]
    classB = classes[1]
    
    Da = D[L==classA,:]
    Db = D[L==classB,:]
    
    rc = stats.ttest_ind(Da, Db, axis=0, equal_var=False)
    ttest_scores = sp.absolute(rc[0])
    ttest_pvals = rc[1]
    TT = sp.vstack((sp.array(ttest_scores),sp.array(ttest_pvals))).T
    
    y = zip(ttest_scores,ttest_pvals,range(len(ttest_scores)))
    y = sorted(y,reverse=True)
    
    if return_mat:
        return y, TT
    else:
        return y

if __name__ == '__main__':
    pass