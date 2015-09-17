'''
Created on Jun 3, 2013
@author: Stephen O'Hara

Functions in support of building classifiers on genomic data
using a combination of typical data-driven machine-learning
techniques combined with information from ontologies, such as
the Gene Ontology (GO).
'''
import ifr
import scipy as sp
from ifr.common.constants import IFR_INFLUENZA_FILE_GENES, IFR_INFLUENZA_FILE_ACC, IFR_LUNG_FILE_GENES, IFR_LUNG_FILE_ACC
from ifr.common.constants import flu_genelists

PC_KNOWN_ALIASES = {'SMCX':'JARID1C', 'CRK7':'CRKRS', 'SMCY':'JARID1D', 'M26880':'UBC', 'NOD2':'IBD1'}

#Gene names that can't be resolved back to an index in our features.
PC_KNOWN_BAD = ['HSAJ2425']

class Pathway_Classifier():
    '''
    Instances of this class are used to build classifiers based upon
    features selected with a combination of IFR and ontology analysis
    '''
    def __init__(self, traindat, testdat, feature_names, numIter=30, ifr_obj=None):
        '''
        Constructor.
        @param traindat: Training data specified as a tuple (D,L) where D is the data matrix, rows are samples
        and columns are features. L is a label vector.
        @param testdat: Tuple (D,L), as per traindat, but this specifies the testing
        data.
        @param feature_names: A list which has a gene name for each feature column. Length of feature_names
        should be same as number of columns in the training data matrix.
        @param ifr_obj: If an iterative feature removal object already exists for this data, you may
        provide it here. Otherwise, specify None, and a new ifr object can be built internally.
        '''
        self.traindat = traindat
        self.testdat = testdat
        self.feature_names = feature_names
        self.numIter = numIter
        self.ifr = ifr_obj
        self.ifra = None
        self.pathways = None
        self.pathway_accuracies = None
        self.pathway_pair_accuracies = None
                
        #a dict to cache known GATHER gene name -> name as indexed in gene list
        #We start with one or more known-aliases that are especially problematic to automatically
        #resolve with web queries. Possibly because the Duke Data has a few obsolete gene labels.
        self.known_aliases = PC_KNOWN_ALIASES
        self.known_bad = PC_KNOWN_BAD
        
    def initAnalysis(self, IFR_Class=ifr.IFR, IFRA_Class=ifr.IFR_Analysis_GO, IFR_kwargs={}, IFRA_kwargs={} ):
        '''
        Use this method to specify the IFR class and parameters as well as the IFR Analysis class
        and parameters used to construct the pathway classifiers. 
        '''
        if self.ifr is None:
            print "===================================="
            print "Performing Iterative Feature Removal"
            print "===================================="
            self._buildIFR(IFR_Class, **IFR_kwargs)
            
            
        if self.ifra is None:
            print ""
            print "===================================="
            print "Performing Ontology Analysis of IFR"
            print "===================================="
            self._buildIFRA(IFRA_Class, **IFRA_kwargs)
        
        
    def _buildIFR(self, IFR_Class = ifr.IFR, **ifr_kwargs ):
        '''
        Builds an iterative feature removal object to be used to select the features
        which are then used to find discriminative pathways.
        @param IFR_Class: Defaults to the base IFR class, but subclasses could work as well
        @param ifr_kwargs: A dictionary of keyword arguments that will be passed onto
        the constructor of the ifr object. If None, then a default set is used that
        is appropriate for the influenza data.
        '''
        (D,T) = self.traindat
        (D2,T2) = self.testdat 
        if len(ifr_kwargs) < 1:
            ifr_kwargs = {'C':1.0, 'common_scale':False}
            
            #def svm_engine(Train, Test, perm=None, common_scale=True, verbose=False, no_normalization=False,
            #   loss='l2', penalty='l1', dual=False, C=0.5, tol=0.0001, **kwargs):
            
        ifrx = IFR_Class( (D, T), (D2,T2), feature_names=self.feature_names, **ifr_kwargs)       
        ifrx.computeAllIterations()
        self.ifr = ifrx
        
    def _buildIFRA(self, IFRA_Class = ifr.IFR_Analysis_GO,  **kwargs):
        '''
        An IFRA object is an IFR Analysis object, which encapsulates the querying of
        an external ontology (GO, KEGG) given the features selected at each IFR iteration.
        @param IFRA_Class: The name of the class for the IFR Analysis, the default being
        ifr.IFR_Analysis_GO
        @param kwargs: The keyword arguments will be passed onto the compute_annotations_per_iteration()
        method of the IFRA_Class instance.
        '''
        if self.ifr is None:
            print "You must have an IFR object before performing IFR Analysis."
            print "You may need to call the _buildIFR() method prior to this method."
            return
        
        ifra = IFRA_Class(self.ifr)
        ifra.compute_annotations_per_iteration(num_iterations=self.numIter, **kwargs)
        self.ifra = ifra
            
    def _pathway_classification(self,genelist,common_scale=False):
        '''
        Given a genelist from the data set, such as a collection of genes belonging to a specific
        pathway, this function computes the accuracy of a standard (L2) SVM classifier
        '''
        (D,L) = self.traindat
        (D2,L2) = self.testdat
            
        G = self.feature_names
        idxs = []
        for x in genelist:
            while x in self.known_aliases: x = self.known_aliases[x]
            #alias could be key to another alias...so go down until we stop finding substitutions
            
            if x in self.known_bad: continue #skip the ones we know are not resolvable
            
            try:
                idx = G.index(x)
            except(ValueError):                 
                print "Gene %s is not a known feature name. Querying for official name..."%x
                xn = ifr.get_official_name(x, entrezgene=True)
                if (xn is None) or (xn==x):
                    xn = None
                    self.known_bad.append(x)
                    print "No other official name found, gene will be skipped."
                else:                    
                    print "Substituting %s for %s."%(xn,x)
                    self.known_aliases[x] = xn
                
                try:
                    idx = G.index(xn) if (not xn is None) else None
                except(ValueError):
                    idx = None  #happens when official name is found, but still not in list
                
            if not idx is None: idxs.append(idx)
        
        rc = ifr.svm_engine((D,L), (D2,L2), perm=idxs, common_scale=common_scale, verbose=False,
                            no_normalization=False, loss="L2", penalty="L2", C=1.0)
        return rc
        
    def getGenelistForPathway(self, pw_name):
        '''
        Use this function to get the list of genes that are part of a given
        pathway classifier. You must have already computed pathway classifiers.
        This function should be used instead of directly accessing the .pathways
        member variable because this will replace the known aliases.
        '''
        gl1 = self.pathways[pw_name]
        gl2 = [ (self.known_aliases[g] if g in self.known_aliases else g) for g in gl1]
        return gl2
        
        
    def computePathwayClassifiers(self, numIter=None, min_size=1):
        '''
        Determines which pathways are represented by the genes selected within numIter iterations
        of Iterative Feature Removal, and then builds a classifier for each pathway to determine
        how discriminative each pathway is.
        @param numIter: How many iterations of the IFR gene sets should be used in constructing the
        'pathways' used to build the pathway classifiers. This number should be <= self.numIter,
        which is the number of iterations that annotations were processed for. If None, then
        all iterations with annotations will be processed.
        @param min_size: What is the minimum number of genes in an IFR iteration that must
        have the pathway label for that pathway to be used.
        @return (acc_dict, d), where acc_dict is a dictionary keyed by pathway name and value
        is the classification accuracy; d is a dictionary keyed by pathway name with values being
        the list of genes that were included in that pathway from the GO analysis.
        '''
        assert not self.ifra is None
        if numIter is None:
            #all iterations in the IFRA are used
            numIter = self.numIter
        
        d = self.ifra.annotationMembership(numIter=numIter)
        d2 = { k:v for (k,v) in d.items() if len(v) >= min_size}
        print "There are %d annotations with at least %d members."%(len(d2), min_size)
        
        acc_dict = {}

        for (k,gl) in d2.items():
            rc = self._pathway_classification(gl)
            acc_dict[k] = rc[0]
            print "%s %3.2f"%(k, rc[0])
            
        self.pathways = d2
        self.pathway_accuracies = acc_dict
        return (self.pathways, self.pathway_accuracies)
    
    def computePathwayPairsClassifiers(self, verbose=False):
        '''
        After having computed the pathway classifiers, this method will
        compute the accuracy when combining the genes from pairs of pathways.
        @return: A results dictionary keyed by (i,j), which is the combination
        of pathways indexed by i and j, and the values are the tuple:
        (accuracy, pathway_name_i + pathway_name_j, number_of_genes)
        '''
        assert not self.pathways is None
        
        N = len(self.pathways)
        keys = self.pathways.keys()
        acc = {}
        
        for i in range(N):
            for j in range(i+1,N):
                p1 = keys[i]
                p2 = keys[j]
                if verbose: print (i,j), "%s + %s"%(p1,p2)
                gl1 = self.pathways[p1]
                gl2 = self.pathways[p2]
                gl = gl1 + gl2
                rc=self._pathway_classification(gl)
                
                #store the results
                acc[(i,j)] = (rc[0], "%s + %s"%(p1,p2), len(gl1)+len(gl2) )
            
        self.pathway_pair_accuracies = acc
        return self.pathway_pair_accuracies
   
    def getTopPathways(self, X=10):
        '''
        Returns the top X highest-accuracy pathways
        @param X: How many of the top discriminative pathways do you want returned. Specify None for
        all.
        @return: (pathway_name, accuracy) sorted by accuracy (high to low)
        '''
        assert not self.pathway_accuracies is None
        tmp = sorted( [(v,k) for (k,v) in self.pathway_accuracies.items()], reverse=True)
        if X is None:
            return tmp
        else:
            return tmp[0:X]
        
    def getTopPathwayPairs(self, X=10):
        assert not self.pathway_pair_accuracies is None
        a = self.pathway_pair_accuracies.values()
        tmp = sorted(a, reverse=True )
        if X is None:
            return tmp
        else:
            return tmp[0:X]
        

    def getPathwayPairsDict(self):
        '''
        Return a dictionary where the key is a pathway pair string,
        and the value is the accuracy.
        Example: { 'cell cycle + protein metabolism':0.68 }
        '''
        assert not self.pathway_pair_accuracies is None
        tmp = self.pathway_pair_accuracies.values()
        d = { pwp:acc for (acc, pwp, _) in tmp }
        return d
        

def test_pathway_classifier(numIter=10):
    (D,L,_,_) = ifr.load_flu_mat()
    (D2,L2,_,_) = ifr.load_H1N1_mat()
    G = ifr.load_gene_ids(short_name_only=True)
    pc = ifr.Pathway_Classifier((D,L), (D2,L2), feature_names=G, numIter=numIter)
    pc.initAnalysis()
    _rc=pc.computePathwayClassifiers(numIter=numIter)
    _rc2=pc.computePathwayPairsClassifiers()
    print "Top 10 pathway pairs:"
    print pc.getTopPathwayPairs(10)
    
    return pc

    
def pathway_classification(genelist, feature_names=None, traindat=None, testdat=None, return_err_idxs=False, 
                           common_scale=False, verbose=True ):
    '''
    Given a genelist from the flu data set, such as a collection of genes belonging to a specific
    pathway, this function computes the accuracy of a standard (L2) SVM classifier
    @param genelist: The list of genes forming the pathway
    @param feature_names: The list of feature (gene) names associated with the columns of test/train data.
    The feature_names must contain the genes in genelist, or a lookup will be performed to find the matching alias.
    If None, then feature names will be the genes from the Duke Influenza data.
    @param traindat: If None, then H3N2 data will be used. Else, specify the tuple (D,L) where D is the
    data matrix (samples in rows) and L is the label vector
    @param testdat: If None, then H1N1 data will be used. Else, specify (D2,L2) tuple for test data.
    @param return_err_idxs: If true, then the indexes in the test set where the classifier is wrong will
    be returned.
    @param common_scale: If true, then the test data will be scaled using the training data mean/std, else
    it will be scaled using its own mean/std.
    @param verbose: If true, more output will be displayed.
    @return: Either returns rc or rc, err_set, where rc is the return from the svm engine, which is the tuple
    (test_accuracy, factors, clf, train_accuracy). See L{ifr.svm_engine}.
    '''
    global PC_KNOWN_BAD
    global PC_KNOWN_ALIASES
    
    if traindat is None:
        (D,L,_,_) = ifr.load_flu_mat()
    else:
        (D,L) = traindat
        
    if testdat is None:
        (D2,L2,_,_) = ifr.load_H1N1_mat()
    else:
        (D2,L2) = testdat
        
    G = feature_names if (not feature_names is None) else ( ifr.load_gene_ids(short_name_only=True) )
    idxs = []

    for x in genelist:
        while x in PC_KNOWN_ALIASES:
            if verbose: print "%s is known alias to %s."%(x, PC_KNOWN_ALIASES[x])
            x = PC_KNOWN_ALIASES[x]
            #alias could be key to another alias...so go down until we stop finding substitutions
        
        if x in PC_KNOWN_BAD:
            if verbose: print "%s is known bad, skipping gene."%x
            continue       
        
        try:
            idx = G.index(x)
        except(ValueError): 
            #x is probably an alias to a gene name in G
            if verbose:
                print "Gene %s is not a known feature name. Querying for official name..."%x
            xn = ifr.get_official_name(x, entrezgene=True)
            if (xn is None) or (xn==x):
                if verbose: print "No other official name found, gene will be skipped."
                xn = None
                PC_KNOWN_BAD.append(x)
            else:                    
                print "Substituting %s for %s."%(xn,x)
                PC_KNOWN_ALIASES[x] = xn
             
            try:
                idx = G.index(xn) if (not xn is None) else None
            except(ValueError):
                idx = None  #happens when official name is found, but still not in list
            
        if not idx is None: idxs.append(idx)

    rc = ifr.svm_engine((D,L), (D2,L2), perm=idxs, common_scale=common_scale, verbose=False,
                        no_normalization=False, loss="L2", penalty="L2", C=1.0)
    
    if return_err_idxs:
        IX = sp.array( range(len(L2)))
        errs = IX[ rc[4] != L2]  #indexes where prediction is not correct
        err_set = set(errs)
        return rc, err_set
    else:
        return rc

def all_GO_annotations_classification(ifra=None, depth=5, numIter=30):
    '''
    Compute the classification accuracies for the set of genes corresponding to
    GO annotations at a given depth in the ontology.
    @param ifra: An instance of IFR_Analysis_GO object that has already computed the annotations
    per iteration at a given ontology depth. Set to None to have this function create a new one.
    @param depth: This parameter is ignored if ifra is not none, in which case ifra.depth is the
    ontology depth. If ifra is none, then this is the ontology depth used when creating a new
    analysis object.
    @param numIter: How many iterations of the IFR gene sets should be used in constructing the
    'pathways' used to build the pathway classifiers.
    @return: (acc_dict, d), where acc_dict is a dictionary keyed by pathway name and value
    is the classification accuracy; d is a dictionary keyed by pathway name with values being
    the list of genes that were included in that pathway from the GO analysis.
    '''
    if ifra is None:
        ifra = ifr.IFR_Analysis_GO(IFR_INFLUENZA_FILE_GENES, ontology_depth=depth)
        ifra.compute_annotations_per_iteration(saveas=None, bf_thresh=0.0, num_iterations=numIter, homologs=False)
    
    d = ifra.annotationMembership(numIter=numIter)
    d2 = { k:v for (k,v) in d.items() if len(v) >= 5}
    print "There are %d annotations with at least 5 members."%len(d2)
    
    acc_dict = {}
    
    (D,L,_,_) = ifr.load_flu_mat()
    (D2,L2,_,_) = ifr.load_H1N1_mat()
    for (k,gl) in d2.items():
        rc = pathway_classification(gl, traindat=(D,L), testdat=(D2,L2))
        acc_dict[k] = rc[0]
        print "%s %3.2f"%(k, rc[0])
        
    return acc_dict, d2
    
def all_GO_annotation_pairs_classification(d, verbose=True):
    '''
    Using the results of the all_GO_annotations_classification, this function
    computes the accuracies of all pairs of the pathway classifiers.
    @param d: The dictionary { pathway_name:gene_list, ... }
    @return: A results dictionary keyed by (i,j), which is the combination
    of pathways indexed by i and j, and the values are the tuple:
    (accuracy, pathway_name_i + pathway_name_j, number_of_genes)
    '''
    (D,L,_,_) = ifr.load_flu_mat()
    (D2,L2,_,_) = ifr.load_H1N1_mat()
    
    N = len(d)
    keys = d.keys()
    acc = {}
    
    for i in range(N):
        for j in range(i+1,N):
            p1 = keys[i]
            p2 = keys[j]
            if verbose: print (i,j), "%s + %s"%(p1,p2)
            gl1 = d[p1]
            gl2 = d[p2]
            gl = gl1 + gl2
            rc=pathway_classification(gl, (D,L), (D2,L2))
            
            #store the results
            acc[(i,j)] = (rc[0], "%s + %s"%(p1,p2), len(gl1)+len(gl2) )
        
    return acc
    
def all_pathway_classification():
    '''
    Compute classification rates for all pathways in flu_genelist
    '''
    (D,L,_,_) = ifr.load_flu_mat()
    (D2,L2,_,_) = ifr.load_H1N1_mat()
    G = ifr.load_gene_ids(short_name_only=True)
        
    accListF = []
    accListS = []
    errors = {}
    #for each pathway, get the gene list and compute classification power
    for p in sorted( flu_genelists.PATHWAYS.keys()):
        (_pathway_title, genelist) = flu_genelists.PATHWAYS[p]
        (rc, err_set) = pathway_classification(genelist, feature_names=G, traindat=(D,L), testdat=(D2,L2),
                                               return_err_idxs=True)
        accListF.append( rc[0] )
        accListS.append( "%3.1f"%(100.0*rc[0]))
        errors[p]=err_set
        #print pathway_title, rc[0]
        
    return accListF, accListS, errors
    
def all_pathway_pairs_classification():
    '''
    computes the classification accuracy using all
    pairs of pathways from the flu_genelist SSVM IFR
    '''
    acc = {}  #the classifier accuracy combining pathway i with j
    (D,L,_,_) = ifr.load_flu_mat()
    (D2,L2,_,_) = ifr.load_H1N1_mat()
    G = ifr.load_gene_ids(short_name_only=True)
    pw_dict = flu_genelists.PATHWAYS
    
    #get the mistakes that each individual pathway classifier makes
    (accList,_,errors) = all_pathway_classification()
        
    for i in range(1,13):
        for j in range((i+1),13):
            #find error overlap
            e1 = errors[i]
            e2 = errors[j]
            f1 = len( e1.intersection(e2) )
            #f2 = len( e1.union(e2) )
            overlap_frac = float(f1)/57
            #find accuracy of the two single pathway classifiers
            a1 = accList[i-1]
            a2 = accList[j-1]
            max_acc = max(a1,a2)
            #compute combined accuracy
            (p1,gl1)=pw_dict[i]
            (p2,gl2)=pw_dict[j]
            gl = gl1+gl2
            rc=pathway_classification(gl, feature_names=G, traindat=(D,L), testdat=(D2,L2))
            #(genelist, feature_names=None, traindat=None, testdat=None, return_err_idxs=False, 
            #               common_scale=False, verbose=True ):
            
            #store the results
            acc[(i,j)] = (rc[0], "%s + %s"%(p1,p2), overlap_frac, max_acc )
                                      
    return acc    

def all_IFR_pairs_classification(ifr, numIter=20):
    '''
    computes the classification accuracy using all
    pairs of IFR iterations from the flu_genelist SSVM IFR
    '''
    acc = {}  #the classifier accuracy combining iteration i with j
    (D,L,_,_) = ifr.load_flu_mat()
    (D2,L2,_,_) = ifr.load_H1N1_mat()
    
    removed = ifr.get_removed_features()
    
    #compute L2 SVM test accuracies for all iterations upto numIter
    print "Computing L2 SVM test accuracies for each IFR iteration."
    test_acc_list = []
    for x in range(numIter):
        glx = removed[x]
        tmp=pathway_classification(glx, (D,L), (D2,L2))
        test_acc_list.append(tmp[0])
    
    print "Computing L2 SVM test accuracies for all pairs of IFR iterations."
    cur = 0
    total = (numIter/2)*(numIter-1)    
    for i in range(numIter):
        for j in range((i+1),numIter):
            ifr.print_progress(cur, total)
            cur+=1
            
            gl1 = removed[i]
            gl2 = removed[j]
            gl = gl1+gl2
            
            a1 = test_acc_list[i]
            a2 = test_acc_list[j]
            max_acc = max(a1,a2)
            
            #compute combined accuracy
            rc=pathway_classification(gl, (D,L), (D2,L2))
            
            #store the results
            acc[(i,j)] = (rc[0], "IFR %d + IFR %d"%(i,j), a1, a2, max_acc )
    
    res = sorted( acc.values(), reverse=True)
                              
    return res, test_acc_list