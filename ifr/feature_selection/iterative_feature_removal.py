'''
Created on Dec 17, 2012
@author: Stephen O'Hara

Feature selection via iterative feature removal.
Each iteration trains a sparse model which selects
only a small number of features to fit the data.
The test accuracy is noted, and then the selected
features are removed from the data set, and the next
iteration occurs. The subsequent iteration's model
must fit the data with a new set of features...

The idea is to determine when a good model can no
longer be fit to the data, and use those features
removed at all iterations prior to this step.
'''

import ifr
import scipy as sp
import cPickle
import pylab as pl
import sys
import math

IFR_PARTITION_RAND = 'Random'
IFR_PARTITION_SUBJ = 'Subject'

#Visualization types for plot method
IFR_PLOT_HEATMAP = 0
IFR_PLOT_BARCHARTS = 1
    
class IFR:
    '''
    IFR stands for Iterative Feature Removal,
    and is a mechanism for selecting the maximum-sized
    relevant subset of features in a data set using
    successive iterations of removing the best K features
    as selected by a sparse machine learning method.
    '''    
    def __init__(self, traindat, testdat, feature_names=None,
                 engine=ifr.ML_ENGINE_SVM, verbose=False, **ml_kwargs):
        '''
        constructor
        @param traindat: The training data over which feature selection
        is performed. Specify as a tuple (D,T) where D is the data matrix,
        samples are in rows, and T is a vector of integer labels. This current
        version only supports binary classification, so there should be only
        two unique label values in T.
        @param testdat: The testing (or validation) data, which is not used
        in feature selection, but is used to measure the accuracy of the selected
        features on the withheld data. Same format as traindat, a tuple (D,T)
        @param feature_names: Optional. If provided, this is a list of strings where each
        element of the list is the name for the feature of the corresponding column
        @param engine: The sparse machine learning engine used
        to train models and select non-zero coefficients. See the
        ml_engines.py module in the "common" package for the interface.        
        @param ml_kwargs: Keyword arguments to be passed to machine learning engine, such
        as parameter settings, stopping criteria, etc.
        @note: partition type and test frac parameters are only used if
        no fixed testing data set was provided (testdat is None). If during object construction,
        a fixed test set was specified, then these parameters are ignored.
        '''
        (D,T) = traindat
        self.D = D
        self.N = D.shape[0]  #number of samples
        self.P = D.shape[1]  #number of columns/features
        self.T = T 
        
        (D2,T2) = testdat
        self.Dt = D2
        self.Tt = T2
        
        self.feature_names = feature_names
        self.verbose = verbose
        self.ml_kwargs = ml_kwargs
        self.engine = engine
        self.iteration = 0
        self.removed = [] #list of lists. Each list consists of col idxs selected during that iteration.
        self.removed_count = 0
        self.train_acc_list = []
        self.test_acc_list = []
        
    def __getitem__(self, n):
        if self.feature_names is None:
            return self.removed[n]
        else:
            idxs = self.removed[n]
            if type(idxs[0]) == list:
                #n was a slice, not a single index, so a list of lists is returned
                rc = [ [(x,self.feature_names[x]) for x in idxList] for idxList in idxs]
            else:
                #n is a single index, so a single list is returned
                rc = [ (x, self.feature_names[x]) for x in idxs]
            return rc
        
    def __iter__(self):
        return self
        
    def next(self):
        self._iterate()
        return self.removed[-1]
        
    def computeAllIterations(self):
        '''
        Will perform iterative removals until there are no more features to be removed.
        '''
        #this object is an interator, and this will stop when there are no more iterations to process
        cur_count = 1
        for _ in self:
            print ".", 
            if cur_count % 40 == 0: print ""       	
            sys.stdout.flush()
            cur_count +=1
        
    def getRollingAvgTestAcc(self, window_size=10):
        return ifr.rolling_average(self.test_acc_list, window_size=window_size)
    
    def getRollingAvgTrainAcc(self, window_size=10):
        return ifr.rolling_average(self.train_acc_list, window_size=window_size)
    
    def getNumIterations(self):
        return len(self.removed)
    
    def getRemovedCount(self, numIter=None):
        '''
        Returns the total number of genes removed over the iterations which
        have currently been processed.
        @param numIter: If None, this function simply returns the value
        of the member variable self.removed_count. Else, if a number N is
        provided, this method returns the count of genes removed for the
        first N iterations.
        '''
        if numIter is None: return self.removed_count
        
        assert len(self.removed) >= numIter
        
        ct = 0
        for r in self.removed[0:numIter]:
            ct += len(r)
        
        return ct
    
    def getRemainingCount(self):
        return self.P - self.removed_count
    
    def getAllRemovedFeatures(self):
        '''
        @return: A flattened list of all col indxs in self.removed
        '''
        all_removed = []
        for lst in self.removed:
            all_removed += lst            
        return all_removed
    
    def getAllRemainingFeatures(self):
        '''
        @return: a list of all col idxs that have NOT been removed yet
        '''
        set1 = set( range(self.P) )
        set2 = set( self.getAllRemovedFeatures() )
        return sorted( list( set1 - set2 ) )
    
    def save(self, fn="IFR_results.p"):
        '''
        Saves the entire ifr data object to a pickle file
        '''
        cPickle.dump(self, open(fn,"wb"))
        print "Saved IFR object to %s"%fn
    
    def get_removed_features(self):
        '''
        Returns the iterative removals using the
        feature names instead of the indexes, as is
        stored in self.removed
        '''
        assert not self.feature_names is None
        results = []
        for itr in self.removed:
            featurelist = [ self.feature_names[idx] for idx in itr]
            results.append(featurelist)
        return results          
        
        
    def export_iterations(self, fn="iterations.txt"):
        '''
        Saves a text file which has the contents of self.removed,
        where each row is an iteration, and has the common names
        of the features separated by commas in each row.
        '''
        with open(fn,"w") as f:
            for itr in self.removed:
                if not self.feature_names is None:
                    featurelist = [ self.feature_names[idx] for idx in itr]
                else:
                    featurelist = [ str(idx) for idx in itr ]
                lineout = ",".join(featurelist)
                f.write("%s\n"%lineout)
        print "Saved iterative feature removal data to: %s"%fn               
    
    def plotResults(self, titlestr="", ylimits=[0.5,1.05], plotfunc = pl.semilogx, ylimitsB=[0,101],
                     legend_loc=3, show=True ):
        pl.figure(num=None, figsize=(15,5))
        xvals = range(1, (1+len(self.removed)) )
        #Two subplots. One the left is the test accuracy vs. iteration
        pl.subplot(1,2,1)
        plotfunc(xvals, self.test_acc_list, "b", label="Test Accuracy")    
        pl.hold(True)
        plotfunc(xvals, self.getRollingAvgTestAcc(window_size=10), "r", label="Test Acc (rolling avg)")
        plotfunc(xvals, self.getRollingAvgTrainAcc(window_size=10), "g--", label="Train Acc (rolling avg)")
        pl.ylim(ylimits)
        if titlestr == "":
            pl.title("Iterative Feature Removal")
        else:
            pl.title(titlestr)
        pl.ylabel("Test Accuracy")
        pl.xlabel("Iteration")
        pl.legend(loc=legend_loc) #3=lower left
        pl.hold(False)
        
        #second subplot. On the right is the number of features removed per iteration
        pl.subplot(1,2,2)
        Ns = [ len(lst) for lst in self.removed ]
        pl.semilogx(xvals, Ns, "bo", label="#Features per Iteration")
        pl.xlabel("Iteration")
        pl.ylabel("Number of Features Selected")
        pl.title("Number of Features Removed per Iteration")
        pl.ylim(ylimitsB)
        
        pl.subplots_adjust(left=0.05, bottom=0.15, right=0.95, top=0.90, wspace=0.20, hspace=0.20)
        if show: pl.show()
        
    def plotSolo(self, titlestr="", ylimits=[0.25,1.05], plotfunc = pl.semilogx, legend_loc=3, show=True ):
        '''
        Much the same as plotResults() method, but only has the 'main' IFR accuracy curve. This one
        is better if you want to embed several IFR curves in a single figure.
        '''
        xvals = range(1, (1+len(self.removed)) )
        
        plotfunc(xvals, self.test_acc_list, "b", label="Test Accuracy")    
        pl.hold(True)
        plotfunc(xvals, self.getRollingAvgTestAcc(window_size=10), "r", label="Test Acc (rolling avg)")
        plotfunc(xvals, self.getRollingAvgTrainAcc(window_size=10), "g--", label="Train Acc (rolling avg)")
        pl.ylim(ylimits)
        if titlestr == "":
            pl.title("Iterative Feature Removal")
        else:
            pl.title(titlestr)
        pl.ylabel("Test Accuracy")
        pl.xlabel("Iteration")
        if not legend_loc is None: pl.legend(loc=legend_loc) #3=lower left
        pl.hold(False)
         
        if show: pl.show()
        
    def _iterate(self):
        '''
        Perform one iterative removal. Note, the intended way to use this
        class after construction is to call the next() method in a loop.
        This is an internal method called by next().        
        '''
        cur_features = self.getAllRemainingFeatures()
        if len(cur_features) < 1:
            if self.verbose: print "All features have been removed."
            raise StopIteration
        else:
            if self.verbose: print "Num available features:", len(cur_features)
        
        #Data filtered for only remaining columns
        D2 = self.D[:,cur_features]
        trn = (D2, self.T)
        tst = (self.Dt[:,cur_features], self.Tt)
        
        (test_acc, features_used_in_model, _, train_acc, _) = self.engine( trn, tst, **self.ml_kwargs )
        if self.verbose:
            print "Train accuracy: %f"%train_acc
            print "Test accuracy: %f"%test_acc
            print "Num features: %d"%len(features_used_in_model)
        
        #remember, there is an inverse mapping (lookup) required to get the feature column from
        # self.D (unfiltered list) from the non-zero coeff indxs learned from D2
        selected_features = self._getFeatureIdxs(cur_features, features_used_in_model)
        if self.verbose:
            print "Selected features in this iteration:"
            print selected_features            
            
        if len(selected_features) < 1:
            if self.verbose: print "No selected features in remaining data. Stopping."
            raise StopIteration
        
        self.removed.append(selected_features)
        self.removed_count += len(selected_features)
        self.train_acc_list.append(train_acc)
        self.test_acc_list.append(test_acc)
        
    def _getFeatureIdxs(self, features_set, model_idxs):
        '''
        The model returns a list of columns with non-zero coefficients. However, the model was
        trained on an iteration where some of the columns from the original data were removed.
        So feature indexes returned by the model are not aligned with the columns of self.D,
        instead we need to determine the columns by referencing the input feature set used to
        filter the data.
        '''
        idxs = [ features_set[i] for i in model_idxs]
        return idxs
        
   
def test_ifr(max_iter=600, plot=True):
    '''
    test ifr iterations using a fixed test data set
    H3N2 data where half of all subjects (not rows) are partitioned
    into the test set.
    @param max_iter: The maximum number of iterations to process.
    '''
    (D,T,_,S) = ifr.load_flu_mat()
    (Dtrain, Ttrain, Dtest, Ttest, test_ids) = ifr.subject_partition_proportional(D, T, S, frac=0.50)
    print "Test ids: %s, with labels: %s"%(str(test_ids),str(Ttest))
    genes = ifr.load_gene_ids(short_name_only=True)
    ifrx = IFR( (Dtrain, Ttrain), (Dtest,Ttest), feature_names=genes )
    for i,x in enumerate(ifrx):  #IFR is an iterator, each iteration returns the features at that iteration
        print "Iteration %d: %d Features Selected. %d Remaining."%((i+1),len(x),ifrx.getRemainingCount())
        if i>=max_iter: break
    
    print "Features from first iteration:"
    print ifrx[0]  #note IFR class supports __getitem__, so results for any iteration processed can be returned
                    # including slices like ifr[0:10] for the first 10 iterations, list of lists returned.
    
    if plot: ifrx.plotResults(titlestr="IFR: H3N2 Fixed 50 Pct Subject Partition")
    return ifrx


def ifr_load_removals_from_file(fn, G=None):
    '''
    Loads the iterative removal data from a *.csv text file.
    @param fn: The full file name to be loaded. Should be a text csv file
    with one row per removal iteration. Each row should have genes separated by commas.
    @param G: If None, then the returned list-of-lists will reflect exactly what
    is in the .csv file. Otherwise, G is an indexed list of gene names used
    to transform the .csv data as it is loaded. If the .csv data contains integers,
    then G[idx] will replace idx, for all idx in the .csv file. Else, if the .csv
    contains gene names, then G.index(g) will replace g for all g in the .csv file,
    where g is a gene name (string).
    '''
    with open(fn,"r") as f:
        tmp_removals = f.readlines()
    
    rem1 = [ x.strip().split(",") for x in tmp_removals]
    removals = []
    if G is None:    
        for row in rem1:
            rowdat = [ g.strip() for g in row ]
            removals.append(rowdat)
    else:
        #inspect first element to determine if
        # we need to replace indexes in csv file
        # with associated gene names, or if the
        # csv file has gene names which we want
        # to replace with indexes.
        elem1 = rem1[0][0]
        try:
            int(elem1)
            flag=True
        except:
            flag=False
            
        for row in rem1:
            if flag:
                rowdat = [ G[int(g)] for g in row ]                
            else:
                rowdat = [ G.index(g) for g in row ]
            removals.append(rowdat)
    
    return removals

class IFR_Analysis_GO:
    '''
    Class to support the analysis of IFR iterations
    using the publicly available Gene Ontology (GO) using
    the gather web service interface.
    '''        
    def __init__(self, ifr_source):
        '''
        constructor
        @param ifr_source: Provide either a string, which is interpreted as the
        full name/path of a csv text file specifying iterative removal data, OR
        provide an instance of the IFR class.
        '''
        if isinstance(ifr_source, ifr.IFR):
            self.removed = ifr_source.get_removed_features()
        elif type(ifr_source) == str:
            self.removed = ifr.ifr_load_removals_from_file(ifr_source)
        else:
            raise ValueError("Parameter ifr_source must be either an instance of IFR or a filename string.")
    
        self.annotation_dat = None #initialize member variable
     
    def getAnnotationDat(self):
        '''
        @return: The annotation data which is computed after calling GO_annotations_per_iteration()
        method.
        '''
        return self.annotation_dat
    
    def load_annotations_per_iteration(self, fn):
        '''
        Loads a pickle file containing the annotations per iteration, such as that
        which can be saved when calling compute_annotations_per_iteration. This is a time-saver
        when you have valid annotations computed, but just want to re-draw the graphs.
        @param fn: The pickle file which has a single variable called res stored in it.
        '''
        self.annotation_dat = cPickle.load(open(fn,"rb"))
        print "Pickle file %s has been loaded. Annotations loaded for %d iterations."%(fn,len(self.annotation_dat))
             
    def compute_annotations_per_iteration(self, depth=5, saveas=None, bf_thresh=0.0,
                                      num_iterations=10, homologs=False):
        '''
        Given a file with iterative gene removal data, this function will plot
        a bar chart of the GO annotations at a given level (depth) for each iteration.
        @param depth: The GO ontology depth to use.
        @param saveas: Provide a filename to save the GO annotation results. This will
        be a python pickle file. Specify None to avoid saving the results to a file.
        @param bf_thresh: The bayes factor threshold. An annotation must have at least this
        bayes factor before being counted in the results.
        @param num_iterations: The number of iterations to process
        @param homologs: Select whether to include gene homologs in the analysis.
        '''
        iters = self.removed[:num_iterations]
            
        res = []
        
        print "Starting analysis..."
    
        for i,genelist in enumerate(iters):
            ifr.print_progress(i, len(iters))
            bcdat = self._gen_bar_chart_data(genelist, depth, bf_thresh=bf_thresh, homologs=homologs)
            res.append((i,bcdat,genelist))
        print ""
            
        if not saveas is None:
            cPickle.dump(res, open(saveas,"wb"))
            print "Results saved to file: %s"%saveas
        
        self.annotation_dat = res
        self.depth = depth  #for convenience, remember the depth used to generate the analysis

    def _gen_bar_chart_data(self, genelist, depth, bf_thresh=0.0, homologs=True):
        '''
        Sets up a gather query, and produces a barchart of the annotations
        returned by gather for a specific depth in the ontology.
        '''
        #gather_ = ifr.gather()
        #gather_.query(genelist,homologs=homologs)
        
        rc = self._compute_GO_annotations_at_depth(genelist, depth=depth, bf_thresh=bf_thresh, homologs=homologs)
        
        if rc is None: return []
        else: return rc[0]

    def _compute_GO_annotations_at_depth(self, genelist, depth, bf_thresh=-10.0, homologs=True, gather_obj=None):
        '''
        Computes the data required to produce a bar chart for the genelist where the x-axis are
        GO Annotations of a specified depth, and the y-axis is the count of genes in genelist that have that annotation.
        '''
        if gather_obj is None:
            gather_ = ifr.gather()
            _ = gather_.query(genelist, homologs=homologs)
        else:
            gather_ = gather_obj
            
        rc = gather_.getDataForAnnotation(key='[%d]'%depth, bf_thresh=bf_thresh)
        if rc is None:
            print "Error: No annotations at depth %d for gene list."%depth
            return None
        
        x = rc[1]
        key1 = "02. AnnotationId"
        key2 = "04. Description"
        key3 = "05. Total Genes With Ann"
        key4 = "14. Genes"
        dat = sorted([ (q[key1],q[key2],q[key3],q[key4]) for q in [ gather_.parseAnnotation(r) for r in x] ])
        #dat will be a list of tuples like: [('GO:0009607', 'response to biotic stimulus', 8),...]
        # sorted by annotation id
                    
        return (dat, gather_)
        
    def annotationMembership(self, numIter=None):
        '''
        Parses the annotation_dat "horizontally", to extract the combined set
        of genes for each unique annotation over the iterations.
        '''
        dat = self.annotation_dat if numIter is None else self.annotation_dat[0:numIter]
        
        annot_dict = {}
        for (_ix,d,_) in dat:
            for (_, an, _, agl) in d:
                gl = agl.split(" ")
                if an in annot_dict:
                    annot_dict[an] += gl
                else:
                    annot_dict[an] = gl
        
        #remove duplicates, which can occur if more than one affy_id maps to the same gene id, for example.
        for k in annot_dict.keys():
            lst = annot_dict[k]
            lst2 = list(set(lst))
            annot_dict[k] = lst2
            
        return annot_dict
        
        
    def _plot_GO_bar_chart(self, iteration, title=None, annotation_set=None, xlim=None,
                          suppress_ylab=False, suppress_desc=False, desc_filter=True,
                          color='lavender', show=False):
        '''
        Creates the bar chart visualization of the dat produced by compute_GO_annotations_at_depth()
        function.
        @param iteration: which iteration of the results to show, zero-based index
        @param title: Used in the plot title
        @param annotation_set: Optional. If provided, then the bar chart will show values
        exactly for this list of annotation ids. If dat does not include the annotation id,
        then a zero-count is shown. If dat has an annotation id not in this set, it will not
        be output. If None, then the annotation set is the set of annotations in dat. This
        is useful if you are using this function to plot subfigs and you'd like to have the
        same set of y-labels for each. If provided, format is a list of tuples [ (annotation_id, desc)]
        @param xlim: Optional. If provided, the x-axis will be limited to [start,end] values.
        @param suppress_ylab: If true, only ticks will be rendered on y-axis labels, no text.
        @param suppress_desc: If true, the text descriptions for the GO annotation will not
        be drawn
        @param desc_filter: If suppress_desc is FALSE and this is TRUE, then descriptions will
        be drawn to the bar chart only for those non-zero annotations.
        @param show: If False, the plot won't be shown. This is useful when this function
        is called in a loop to generate subplots in a master figure.
        '''   
        dat = self.annotation_dat[iteration][1]
        depth = self.depth
        
        if not annotation_set is None:
            #we need to alter the labels/counts/desc lists
            # to be keyed off the desired annotation_set
            (descs, labels) = zip(*annotation_set)
            dat_dict =  { l:(d,c) for (l,d,c,_) in dat}
            counts = []
            for lbl in labels:
                tmp = 0 if not lbl in dat_dict else dat_dict[lbl][1]
                counts.append(tmp)
        else:
            labels = [ l for (l,_,_,_) in dat]
            counts = [ c for (_,_,c,_) in dat]
            descs =  [ d for (_,d,_,_) in dat]
       
        maxX = max(counts) if xlim is None else xlim[1]
        
        ylocations = sp.array(range(len(labels)))+0.5
        if not xlim is None: pl.xlim(xlim)
        pl.barh(range(len(labels)), counts, color=color)
        
        if not suppress_desc:
            for i,txt in enumerate(descs):
                if desc_filter:
                    if counts[i] > 0: pl.text(0.01*maxX, ylocations[i], txt, fontsize=10,
                                              ha='left', va='center')
                else:
                    pl.text(0.01*maxX, ylocations[i], txt, fontsize=10, ha='left', va='center')
        
        if suppress_ylab:   
            pl.yticks(ylocations, ['']*len(labels))
        else:
            pl.yticks(ylocations, labels)
            
        pl.ylim(-0.5, ylocations[-1]+0.5)
        if title is None:
            pl.title("Annotations at Depth %d"%depth)
        else:
            pl.title(title)
        pl.gca().get_xaxis().tick_bottom()
        pl.gca().get_yaxis().tick_left()
        pl.subplots_adjust(left=0.15, right=0.95)
        if show: pl.show()

    
    def _common_Y(self):
        '''
        Internal function to generate a combined set of y-axis labels from
        the results data. Used by _plot_engine, _heatmap_engine
        '''
        res = self.annotation_dat
        tmp = []
        all_counts = []
        for (_,dat,_) in res:
            all_counts += [ c for (_,_,c,_) in dat]
            tmp += [ (d,l) for (l,d,_,_) in dat]
        AS = sorted( list(set(tmp)), reverse=True)   #complete Annotation Set
        max_count = int( sorted(all_counts)[-1] )
        return (AS, max_count)
        
    def plot(self, suptitle, numIter=None, visualization=IFR_PLOT_HEATMAP, **kwargs):
        '''
        Plots the results of the GO analysis
        @param suptitle: The main title of the plot
        @param numIter: Default is None, meaning all iterations will be plotted.
        Otherwise specify the number of iterations to include in the plot.
        @param visualization: Which plotting engine to use? Must be one of
        the valid constants defined in this module
        @param kwargs: Additional kwargs to pass to the plotting engine,
        such as a filter count, etc.
        '''
        if visualization == IFR_PLOT_HEATMAP:
            self._heatmap_engine(suptitle, numIter, **kwargs)
        elif visualization == IFR_PLOT_BARCHARTS:
            self._barplot_engine(suptitle, numIter, **kwargs)
        else:
            raise ValueError("Invalid chart type specified.")
    
    def _barplot_engine(self, suptitle, numIter=None, saveas=None, show=True):
        '''
        code to support a barchart visualization of the results, with one
        bar chart per iteration.
        '''
        res = self.annotation_dat if numIter is None else self.annotation_dat[0:numIter]
        
        pl.figure(num=1)
            
        #AS is complete Annotation Set, max_count is the maximum number
        # of genes that appears in any single annotation entry
        (AS, max_count) = self._common_Y()
        
        for (i, dat, genelist) in res:
            if len(dat) < 1: continue
            pl.subplot(2, math.ceil( len(res)/2.0), i+1)
            title = "Iteration %d (%d Genes)"%((i+1), len(genelist))
            flag1=False if i in [0,5] else True
            
            self._plot_GO_bar_chart(i, title, annotation_set=AS, xlim=[0,max_count],
                                  suppress_ylab=flag1, suppress_desc=False, desc_filter=True,
                                  color="yellow", show=False)  
        pl.subplots_adjust(left=0.1, right=0.95, top=0.90, bottom=0.05, wspace=0.1, hspace=0.2)
        if not suptitle is None:
            pl.suptitle("Ontology Annotations per Iteration: %s"%suptitle)
        if not saveas is None:
            pl.savefig(saveas)   
        if show: pl.show()
    
    def _heatmap_engine(self, suptitle, numIter=None, filter_count=None, saveas=None, show=True, cmap="Blues"):
        '''
        An alternative way to visualize the annotations per iterative via
        a heat-map like construct. The rows are the GO annotations and the columns
        are iterations. The color intensity indicates how much a given annotation
        was present in a given iteration
        '''
        res = self.annotation_dat if numIter is None else self.annotation_dat[0:numIter]
        depth = self.depth
        
        pl.figure(num=1,figsize=(20,8))
    
        #AS is complete Annotation Set, max_count is the maximum number
        # of genes that appears in any single annotation entry
        (AS, max_count) = self._common_Y()
        
        #map is a grid Y=annotation X=iteration M(X,Y) = scaled count (c/max c)
        M = sp.zeros( (len(AS), len(res) ) )
        for (col, dat, _) in res:
            if len(dat) < 1: continue
            for (l,d,c,_) in dat:
                row = AS.index((d,l))
                M[row,col] = c #(c*1.0)/max_count
            
        #filter rows / pathways which dont show up often
        if not filter_count is None:
            assert type(filter_count) == int
            vals = sp.sum(M, axis=1)
            M = M[ vals >= filter_count, :]  #only pathways with at least filter count over all iterations
            AS = [ x for i,x in enumerate(AS) if vals[i] >= filter_count ]
            
        pl.imshow(M, interpolation='nearest', aspect=len(res)*1.0/len(AS), origin='lower')
        #for (l,d) in AS:
        #    row = AS.index((d,l))
        #    pl.text(-0.5, row, d[0:40], color="white", verticalalignment='center', fontsize=9)
        
        (descs, _labels) = zip(*AS)
        descs = [ d[0:30] for d in descs]  #truncate long descriptions
        ylocations = sp.array(range(len(descs)))
        pl.yticks(ylocations, descs, fontsize=9, verticalalignment='center')
        
        pl.set_cmap(cmap)
        pl.xticks(range(len(res)), range(1, len(res)+1))
        pl.ylabel("GO Annotation at Depth %d"%depth, fontsize=16)
        pl.xlabel("Iteration", fontsize=16)
        pl.colorbar(ticks=range(1, max_count+1))       
        if not suptitle is None:
            pl.title("Ontology Annotations per Iteration: %s"%suptitle, fontsize=18)
        if not saveas is None:
            pl.savefig(saveas)   
        if show: pl.show()



class IFR_Analysis_KEGG:
    '''
    Class to support the analysis of IFR iterations
    using the publicly available KEGG pathways database.
    '''        
    def __init__(self, ifr_source):
        '''
        constructor
        @param ifr_source: Provide either a string, which is interpreted as the
        full name/path of a csv text file specifying iterative removal data, OR
        provide an instance of the IFR class.
        '''
        if isinstance(ifr_source, ifr.IFR):
            self.removed = ifr_source.get_removed_features()
        elif type(ifr_source) == str:
            self.removed = ifr.ifr_load_removals_from_file(ifr_source)
        else:
            raise ValueError("Parameter ifr_source must be either an instance of IFR or a filename string.")
    
        self.annotation_dat = None #initialize member variable
     
    def getAnnotationDat(self):
        '''
        @return: The annotation data which is computed after calling GO_annotations_per_iteration()
        method.
        '''
        return self.annotation_dat
     
    def load_annotations_per_iteration(self, fn):
        '''
        Loads a pickle file containing the annotations per iteration, such as that
        which can be saved when calling compute_annotations_per_iteration. This is a time-saver
        when you have valid annotations computed, but just want to re-draw the graphs.
        @param fn: The pickle file which has a single variable called res stored in it.
        '''
        self.annotation_dat = cPickle.load(open(fn,"rb"))
        print "Pickle file %s has been loaded. Annotations loaded for %d iterations."%(fn,len(self.annotation_dat))
        
    def compute_annotations_per_iteration(self, saveas=None, affy_names=True, num_iterations=10, species_code="hsa"):
        '''
        Given a file with iterative gene removal data, this function will generate
        the data showing kegg pathways involved in the gene selections of each annotation
        @param saveas: If None, no output will be saved. Otherwise specify a python pickle file to
        store the output.
        @param affy_names: If True, the strings in the csv file are affy probe_ids, and must
        be converted first to Gene Ids. If false, the strings represent Gene Ids already.
        @param num_iterations: The number of iterations to process
        '''
        if affy_names:
            affy_dict = ifr.load_affy_to_geneId_dict()
            
        iters = self.removed[:num_iterations]  
        res = []
        
        print "Starting analysis..."
        
        for i,gl in enumerate(iters):
            #svo_util.print_progress(i, len(iters))
            print "==========================="
            print "Processing iteration %d ..."%(i+1)
            print "==========================="
            genelist = gl if not affy_names else ifr.get_genelist_from_affy(gl, d=affy_dict)
            bcdat = ifr.compute_kegg_annotations(genelist,species_code)
            res.append((i,bcdat,genelist))
        print ""
            
        if not saveas is None:
            cPickle.dump(res, open(saveas,"wb"))
            print "Results saved to file: %s"%saveas
                
        self.annotation_dat = res
    

    def _common_Y(self):
        '''
        Internal function to generate a combined set of y-axis labels from
        the results data. Used by _plot_engine, _heatmap_engine
        '''
        res = self.annotation_dat
        tmp = []
        all_counts = []
        for (_,dat,_) in res:
            all_counts += [ c for (_,_,c) in dat]
            tmp += [ (d,l) for (l,d,_) in dat]
        AS = sorted( list(set(tmp)), reverse=True)   #complete Annotation Set
        max_count = int( sorted(all_counts)[-1] )
        return (AS, max_count)
        
    def plot(self, suptitle, numIter=None, visualization=IFR_PLOT_HEATMAP, **kwargs):
        if visualization == IFR_PLOT_HEATMAP:
            self._heatmap_engine(suptitle, numIter, **kwargs)
        elif visualization == IFR_PLOT_BARCHARTS:
            self._barplot_engine(suptitle, numIter, **kwargs)
        else:
            raise ValueError("Invalid chart type specified.")
    
    def _barplot_engine(self, suptitle, saveas=None, show=True):
        '''
        engine to plot iterations as a series of successive bar charts
        '''
        raise NotImplementedError
        
    def _heatmap_engine(self, suptitle, numIter=None, filter_count=None, saveas=None, show=True, cmap="Reds"):
        '''
        A way to visualize the kegg pathway annotations per iterative via
        a heat-map like construct. The rows are the pathway annotations and the columns
        are iterations. The color intensity indicates how much a given annotation
        was present in a given iteration (number of genes having that annotation)
        '''
        pl.figure(num=1,figsize=(20,8))

        res = self.annotation_dat if numIter is None else self.annotation_dat[0:numIter]
        (AS, max_count) = self._common_Y()
    
        #map is a grid Y=annotation X=iteration M(X,Y) = scaled count (c/max c)
        M = sp.zeros( (len(AS), len(res) ) )
        for (col, dat, _) in res:
            if len(dat) < 1: continue
            for (l,d,c) in dat:
                row = AS.index((d,l))
                M[row,col] = c #(c*1.0)/max_count
            
        #filter rows / pathways which dont show up often
        if not filter_count is None:
            assert type(filter_count) == int
            vals = sp.sum(M, axis=1)
            M = M[ vals >= filter_count, :]  #only pathways with at least filter count over all iterations
            AS = [ x for i,x in enumerate(AS) if vals[i] >= filter_count ]
        
        pl.imshow(M, interpolation='nearest', aspect=len(res)*1.0/len(AS), origin='lower')
        #for (l,d) in AS:
        #    row = AS.index((d,l))
        #    pl.text(-0.5, row, d[0:40], color="white", verticalalignment='center', fontsize=9)
        
        (descs, _labels) = zip(*AS)
        descs = [ d.split(' - ')[0].strip()[0:30] for d in descs]  #truncate long descriptions
        ylocations = sp.array(range(len(descs)))
        pl.yticks(ylocations, descs, fontsize=9, verticalalignment='center')
        
        pl.set_cmap(cmap)
        pl.xticks(range(len(res)), range(1, len(res)+1))
        pl.ylabel("KEGG Pathway", fontsize=16)
        pl.xlabel("Iteration", fontsize=16)
        pl.colorbar(ticks=range(1,max_count+1))     
        if not suptitle is None:
            pl.title("KEGG Pathways per Iteration: %s"%suptitle, fontsize=18)
        if not saveas is None:
            pl.savefig(saveas)   
        if show: pl.show()



if __name__ == '__main__':
    test_ifr()
    
    
    