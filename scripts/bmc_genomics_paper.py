'''
Created on Feb 25, 2013
@author: Stephen O'Hara

Collection of scripts used to produce the graphics and other
derived data in support of our first paper related to our
infectious disease research funding. This paper is provisionally
entitled:

"Iterative Feature Removal Yields Highly Discriminative Pathways",
O'Hara, Wang, Huber, O'Hern, Shattuck, Slayden, Schenkel, and Kirby, submitted 2013
Currently under review by BMC Genomics
'''

import ifr
import os
import scipy as sp
import pylab as pl
import sys
import math
import itertools as it
import cPickle
from scipy.stats import ttest_ind

#from sklearn import svm as svm
    
try:
    import asciitable
except:
    print "Error importing asciitable package. Functions that produce LaTeX tables"
    print " will fail to run."

from ifr.common.constants import IFR_INFLUENZA_FILE_GENES, IFR_INFLUENZA_FILE_ACC, IFR_LUNG_FILE_GENES, IFR_LUNG_FILE_ACC 
from ifr.common.constants import flu_genelists
from ifr.pathways.pathway_classification import all_pathway_classification, all_pathway_pairs_classification


##############################################
#  RESULTS AND DISCUSSION
##############################################

def figure_1_influenza():
    '''
    This function generates the IFR (Iterative Feature Removal) accuracy plots given
    the data produced from running the matlab script which generates the data.
    @param genelist: A csv text file, one row per iteration, each row lists the names of
    the genes with non-zero weights in the model. These are the genes that are removed prior to
    the subsequent iteration.
    @param acclist: A text file in a structured format that provides the test and train accuracy
    for each iteration of the process.
    '''
    
    #THE FOLLOWING CODE WILL COMPUTE THE IFR
    # BUT FOR SPEED, WE SIMPLY LOAD PRE-COMPUTED DATA
    # FROM A SAVED FILE.
    #(D,L,_,_) = ifr.load_flu_mat()
    #(D2,L2,_,_) = ifr.load_H1N1_mat()
    #ifrx = ifr.IFR( (D,L), (D2,L2), common_scale=False)
    #ifrx.computeAllIterations()
    
    #LOAD PRE-COMPUTED IFR DATA FROM A SAVED FILE
    ifrx = _load_stored_IFR('flu')
    
    print "There are %d genes selected in the first 20 iterations."%ifrx.getRemovedCount(numIter=20)
    
    ifrx.plotResults("Influenza", ylimitsB=[0,50])
    return ifrx
    
def figure_1_lung():
    '''
    This function generates the IFR (Iterative Feature Removal) accuracy plots given
    the data produced from running the matlab script which generates the data.
    @param genelist: A csv text file, one row per iteration, each row lists the names of
    the genes with non-zero weights in the model. These are the genes that are removed prior to
    the subsequent iteration.
    @param acclist: A text file in a structured format that provides the test and train accuracy
    for each iteration of the process.
    '''
    
    #PERFORM IFR ON LUNG CANCER DATA
    #(X1,X2) = ifr.load_lung_cancer_dat()
    #ifrx = ifr.IFR( (X1.data, X1.labelsAsIntegers()), 
    #               (X2.data, X2.labelsAsIntegers()), common_scale=True)
    #ifrx.computeAllIterations()
    
    #LOAD PRE-COMPUTED IFR DATA FROM A SAVED FILE
    ifrx = _load_stored_IFR('lung')
    
    print "There are %d genes selected in the first 20 iterations."%ifrx.getRemovedCount(numIter=20)
    
    ifrx.plotResults("Lung Cancer", ylimitsB=[0,40])
    return ifrx

#===========================
# Manual expert analysis of Influenza
#===========================
def figure_2_influenza_pathway_heatmap(numIter=40, cmap="Greens"):
    '''
    Using the curated pathways from Alan, generates
    an iteration-by-iteration heatmap similar to the GO/Kegg figures.
    '''
    assert numIter <= 40  #we only have pathway membership info for genes in first 40 iterations
    
    removals = ifr.ifr_load_removals_from_file(IFR_INFLUENZA_FILE_GENES)
    pw_dict = flu_genelists.PATHWAYS
    
    M = sp.zeros((len(pw_dict), numIter))
    
    descs = []
    for key in pw_dict:
        pw, gl = pw_dict[key]
        descs.append(pw[0:30])
        its = genelist_iterations(gl, removals=removals)
        row = key - 1 #key is 1-based, rows are 0-based
        for col in its:
            if col > numIter: continue  #outside desired range for iterations
            M[row,(col-1)] += 1
    
    #aspect=(numIter*1.0)/len(pw_dict),
    
    pl.imshow(M, interpolation='nearest',  origin='upper')
    ylocations = sp.array(range(len(descs)))
    pl.yticks(ylocations, descs, fontsize=9, verticalalignment='center')
    
    pl.set_cmap(cmap)
    pl.xticks(range(numIter), range(1, numIter+1))
    pl.ylabel("Pathway", fontsize=16)
    pl.xlabel("Iteration", fontsize=16)
    pl.colorbar(ticks=range(1, int(max(M.ravel()))+1))       
    pl.title("Pathway Representation per Iteration (Influenza)", fontsize=18)
    pl.show()            
    return M
    
def table_2_pathways():   
    '''
    Generates the data and latex markup to produce table 2
    Shows the pathways curated by Alan, the genes in each,
    and the classification accuracy of an SVM model built
    using the pathway.
    '''
    data_dict = {}
    ps = []
    ds = []
    #for each pathway, get the gene list and iteration number
    for p in sorted( flu_genelists.PATHWAYS.keys()):
        (pathway_title, genelist) = flu_genelists.PATHWAYS[p]
        gl = sorted(genelist)
        itlist = genelist_iterations(gl)
        tmp = [ "%s (%d)"%(g,i+1) for (g,i) in zip(gl,itlist)]
        outstr = ", ".join(tmp)
        #print pathway_title, outstr
        ps.append("\hline \n %s"%pathway_title)
        ds.append(outstr)
        
    (_,accListS,_) = ifr.all_pathway_classification()
        
    data_dict['Pathway'] = ps
    data_dict['Genes'] = ds
    data_dict['Acc'] = accListS  #classification accuracy of each pathway
    
    cap = r"{\bf Selected pathways relating to genes taken from the first 40 iterations of IFR}"
    cap += " on the Influenza data. The iteration that the genes were discovered is in parentheses."
    cap += " Acc indicates the classification accuracy (percent) using only the genes in each list."
 
    try:
        asciitable.write(data_dict, sys.stdout, Writer = asciitable.Latex,
             col_align="|c|c|p{3.825in}|", latexdict = {'tabletype': 'table*', 
             'preamble': ['\caption{%s}'%cap, r'\label{tbl:pathways}'], #'data_start':r"\hline",
             'data_end':r"\hline", 'header_start':r"\hline" },
             names=['Pathway','Acc','Genes'])
    except:
        print "Error calling asciitable module for output of LaTeX."
    
    return data_dict
    
def table_3_pathway_pairs():
    '''
    script to produce table 3
    Shows that pathways can be combined to
    yield better performing classifiers.
    '''
    acc = ifr.all_pathway_pairs_classification()
    x = sorted( acc.values(), reverse=True )
    top16 = [ ("%2.1f"%(100.0*a), p) for (a,p,_,_) in x[0:16] ]
    (col1, col2) = zip(*top16)
    
    data_dict = {}
    data_dict['Acc'] = col1
    data_dict['Pathway Pair'] = col2
    cap = r"{\bf Classification accuracies improve when combining pathways.} We computed"
    cap += " the classification accuracies when using all pairs of the pathway gene"
    cap += " lists presented in Table 1. Top sixteen results are shown."
    try:
        asciitable.write(data_dict, sys.stdout, Writer = asciitable.Latex,
             col_align="|c|p{4in}|", latexdict = {'tabletype': 'table*', 
             'preamble': [r'\centering','\caption{%s}'%cap, r'\label{tbl:pathway_pairs}'],
             'data_start':r"\hline",
             'data_end':r"\hline", 'header_start':r"\hline" },
             names=['Acc','Pathway Pair'])
    except:
        print "Error calling asciitable module for output of LaTeX."
        
    return data_dict

def figure_3_top17univ_vs_bestclassifier():
    '''
    Generates figure 3 in the manuscript,
    comparison of t-test gene ranking for 17 genes
    selected either by top pathway pair or top
    t-test score.
    '''
    (D,L,_,_) = ifr.load_flu_mat()
    (D2,L2,_,_) = ifr.load_H1N1_mat()
    y = ifr.ttest_features(D, L)
    G = ifr.load_gene_ids(short_name_only=True)
    
    #gene list for best pathway pair classifier    
    gl = flu_genelists.SSVM_BCell + flu_genelists.SSVM_CAM
    #gene list for an equal number of top univariate genes (ttest score)
    gl2 = [ G[idx] for (_,_,idx) in y[0:len(gl)] ]
    
    #rank_gene, list of gene names in ttest ranking order
    rg = [ G[idx] for (_,_,idx) in y ]
    gl_ranks = [ rg.index(g)+1 for g in gl ]  #1-based indexing so best is #1, not #0
    #gl2_ranks = [ rg.index(g)+1 for g in gl2 ]  #this should be 1..len(gl), just to check

    rc1 = ifr.pathway_classification(gl,traindat=(D,L), testdat=(D2,L2))
    rc2 = ifr.pathway_classification(gl2, traindat=(D,L), testdat=(D2,L2))
    
    print "Classification accuracy of selected B-Cell + CAM genes: %3.2f"%rc1[0]
    print "Classification accuracy of equal number of top t-test ranked genes: %3.2f"%rc2[0]
    
    print "T-test ranks of B-Cell + CAM genes:"
    print gl_ranks
    
    pl.subplot(2,1,1)
    gene_expressions_pathway(gl, "B-Cell + CAM", zscore=True, show=False, dat=(D,L))
    pl.subplot(2,1,2)
    gene_expressions_pathway(gl2, "Top 17 T-Test Ranked", zscore=True, show=False, dat=(D,L))
    
    pl.subplots_adjust(left=0.05, bottom=0.15, right=0.95,
                       top=0.95, wspace=0.2, hspace=0.55)
    pl.show()
    
    return (gl, gl2, rc1, rc2, gl_ranks)

#===========================
# Automated analysis of Lung Cancer
#===========================

def figure_4_GO_lung(genelist=IFR_LUNG_FILE_GENES, numIter=30, depth=5, filter_count=3, cmap="Blues"):
    '''
    Generates the heatmap figure showing the GO annotations from depth 5
    at each iteration of IFR over the first 30 iterations
    '''    
    x = ifr.IFR_Analysis_GO(genelist)
    x.compute_annotations_per_iteration(saveas=None, depth=depth, bf_thresh=0.0, num_iterations=numIter, homologs=False)
    x.plot("Lung Cancer", visualization=ifr.IFR_PLOT_HEATMAP, filter_count=filter_count, cmap=cmap)
    return x

def figure_5_top_classifier_generanks():
    '''
    Generates the figure showing the gene ranking distribution
    of the top pathway and pathway pair classifiers compared to
    the optimal SSVM classifier over 4 microarray data sets
    '''
    fn = 'script2_res.p'
    base_dir = os.path.join(ifr.RESULTS_DIR,'discriminative pathways')
    lung_res = os.path.join(base_dir, 'lung', 'lung_%s'%fn)
    flu_res = os.path.join(base_dir, 'influenza', 'flu_%s'%fn)
    prostate_res = os.path.join(base_dir, 'prostate','prostate_%s'%fn)
    bcell_res = os.path.join(base_dir, 'bcell', 'bcell_%s'%fn)
    res_files = (lung_res, flu_res, prostate_res, bcell_res)
    res_list = [ cPickle.load(open(rf,'rb')) for rf in res_files]
    mt_plot_generank_res(res_list, show=True)

def table_6_multidatasets():
    '''
    Generates the values used in the table summarizing x-validation results over
    multiple data sets. This function loads stored results
    to avoid recomputing everything.
    '''
    base_dir = os.path.join(ifr.RESULTS_DIR,'discriminative pathways')
    
    #lung_res = os.path.join(base_dir, 'lung', 'pc_results_lung_80_95.p')
    lung_res = os.path.join(base_dir, 'lung', 'pc_results_lung_80_95_26Oct2013.p')
    
    #flu_res = os.path.join(base_dir, 'influenza','pc_results_flu_80_95.p')
    flu_res = os.path.join(base_dir, 'influenza','pc_results_flu_80_95_26Oct2013.p')
        
    #prostate_res = os.path.join(base_dir, 'prostate','pc_results_prostate_70_75.p')
    prostate_res = os.path.join(base_dir, 'prostate','pc_results_prostate_70_75_26Oct2013.p')
    
    #bcell_res = os.path.join(base_dir, 'bcell','pc_results_bcell_70_75_d5_iter20.p')
    bcell_res = os.path.join(base_dir, 'bcell','pc_results_bcell_70_75_iter20_26Oct2013.p')
    
    
    res_files = [ (lung_res,'Lung'), (flu_res,'Influenza 11-14'), (prostate_res,'Prostate'), 
                 (bcell_res,'BCell Lymph.')]
    ssvm_files = [ os.path.join( os.path.dirname(x), 'ssvm_xval_res.p') for (x,_) in res_files]
        
    print "Loading result files..."
    res_set = [ ( cPickle.load(open(fn,'rb')), ds_name) for (fn,ds_name) in res_files]
    
    print "Loading SSVM cross validation results..."
    ssvm_set= [ cPickle.load(open(fn,'rb')) for fn in ssvm_files]
    
    print "Computing summary results..."
    D = []
    for idx in range( len(res_set)):
        ((pw_res, pwp_res,_), dn) = res_set[idx]
                
        (va,sa,ta,_) = mt_best_validation_score(pw_res, K=1)[0]
        s1 = "%3.1f/%3.1f"%(100*va, 100*sa)
        s2 = "%3.1f"%(100*ta)
        s3 = "%3.1f"%(mt_topK_average(pw_res, K=5)*100)
        
        (va2,sa2,ta2,_) = mt_best_validation_score(pwp_res, K=1)[0]
        s4 = "%3.1f/%3.1f"%(100*va2, 100*sa2)
        s5 = "%3.1f"%(100*ta2)
        s6 = "\\bf{%3.1f}"%(mt_topK_average(pwp_res, K=5)*100)
        
        ssvm_scores = sp.array(ssvm_set[idx]['validation_scores_50trial'])
        top5 = ssvm_set[idx]['test_accuracy_top5']
        va3 = ssvm_scores.mean()
        sa3 = ssvm_scores.std()
        ta3 = top5[0]
        tp5 = sp.mean(top5)
        s7 = "%3.1f/%3.1f"%(100*va3, 100*sa3)
        s8 = "%3.1f"%(100*ta3)
        s9 = "%3.1f"%(tp5*100)
        
        D.append( (dn,s1,s2,s3,s4,s5,s6,s7,s8,s9) )
    D = sp.array(D) #array of strings
    
    data_dict = {}
    col_headers = ('Data Set', 'P Validn', 'P Test', 'P Top 5', 'PP Validn', 'PP Test', 'PP Top 5',
                   'SSVM Vald', 'SSVM Test', 'SSVM Top 5')
    
    for ix, ch in enumerate(col_headers):
        data_dict[ch] = D[:,ix]
    
    cap = r"{\bf Pathway classification on multiple data sets.} More caption here."
    
    asciitable.write(data_dict, sys.stdout, Writer = asciitable.Latex,
             col_align="|l|l|l|l|l|l|l|l|l|l|", latexdict = {'tabletype': 'table*', 
             'preamble': ['\caption{%s}'%cap, r'\label{tbl:multidata_res}'], #'data_start':r"\hline",
             'data_end':r"\hline", 'header_start':r"\hline" },
             names=col_headers)
##############################################
#  Methods
##############################################

def figure_6_SSVM_shelf():
    '''
    SSVM Shelf diagram for methods section
    '''
    fn = os.path.join(ifr.RESULTS_DIR,"IFR_SSVM","h3n2_iter_1_shelf.csv")
    S = sp.loadtxt(fn)
    S = sp.absolute(S)  #magnitudes
    pl.semilogy( range(200), S[0:200], "bo")
    pl.semilogy( range(200), S[0:200], "r--")
    pl.xlabel("Feature Index Sorted by Weight Magnitude")
    pl.ylabel("Weight Magnitude")
    pl.title("200 Largest Weights from an SSVM Model")
    pl.show()


##############################################
#  Supplemental Material
##############################################
    
def figure_8_LR_influenza():
    '''
    IFR with L1 Logistic Regression instead of sparse SVM
    generates the sparse logistic regression
    iterative feature removal graph for use in supplemental section.
    '''
    (D,L,_,_) = ifr.load_flu_mat()
    (D2,L2,_,_) = ifr.load_H1N1_mat()
    G = ifr.load_gene_ids(short_name_only=True)
    ifr_LR = ifr.IFR( (D, L), (D2,L2), engine=ifr.ML_ENGINE_LR,
                      feature_names=G, common_scale=False)
    ifr_LR.computeAllIterations()
    ifr_LR.plotResults("Influenza", ylimitsB=[0,100])
    
def figure_8_LR_lung():
    '''
    generates the sparse logistic regression
    iterative feature removal graph for use in supplemental section.
    '''
    (X1,X2) = ifr.load_lung_cancer_dat()
    D = X1.data; L= X1.labelsAsIntegers()
    D2 = X2.data; L2 = X2.labelsAsIntegers()
    ifr_LR = ifr.IFR( (D, L), (D2,L2), engine=ifr.ML_ENGINE_LR,
                      feature_names=X1.fields)
    ifr_LR.computeAllIterations()
    ifr_LR.plotResults("Lung Cancer", ylimitsB=[0,100])
    
def figure_9_KEGG_lung(genelist=IFR_LUNG_FILE_GENES, numIter=30, filter_count=10, cmap="Reds"):
    '''
    Supplemental figure
    Shows KEGG pathways per iteration on lung cancer data
    '''
    x = ifr.IFR_Analysis_KEGG(genelist)
    x.compute_annotations_per_iteration(affy_names=True, saveas=None, num_iterations=numIter)
    x.plot("Lung Cancer", visualization=ifr.IFR_PLOT_HEATMAP, filter_count=filter_count, cmap=cmap)
    return x
    
def figure_10_univariate_vs_bivariate_separation():
    '''
    Supplemental material
    Two genes, HLA-E and IL17RC, both found in the first 30 iterations of SSVM IFR
    are individually not discriminative in a univariate sense. When combined, however,
    they can partition the training data 90+ % 
    '''
    (D,L,_,_) = ifr.load_flu_mat()
    gl = ['HLA-E', 'IL17RC']
    idxs = [4057, 10635]
    
    Dx = D[:, sp.array(idxs)]
    Dxn, _means, _stds = ifr.normalize_data(Dx)
        
    rc = ifr.svm_engine( (Dxn,L), (Dxn,L), verbose=True, no_normalization=True, loss='l2', penalty='l2', C=1.0)
    
    svm_model = rc[2]
    w = svm_model.coef_[0]
    b = svm_model.intercept_[0]
    print w,b 
    
    xlim=[-2.25,2.25]
    ylim=[-1.75,3.0]
    margin_xs = sp.linspace(xlim[0],xlim[1], 100)
    margin_ys = (-w[0]/w[1])*margin_xs + (-b/w[1])  #equation of line y = mx +b
    
    #plot univariate boxplots of the two genes
    pl.subplot(1,2,1)
    ifr.gene_expression_boxplots(D, L, gl , idxs, zscore=True, titlestr="Univariate Separation of Samples using HLA-E and IL17RC")
    #pl.ylim(ylim)
    
    #plot the bivariate analysis with SVM linear separation
    pl.subplot(1,2,2)
    pl.plot(Dxn[L==1,0], Dxn[L==1,1], "ro")
    pl.plot(Dxn[L==1,0], [ylim[0]+0.05]*27, "rv")  #projection class1 to x-axis
    pl.plot([xlim[0]+0.03]*27, Dxn[L==1,1], "r<")  #projection class1 to y-axis
    
    pl.plot(Dxn[L==-1,0], Dxn[L==-1,1], "bo")
    pl.plot(Dxn[L==-1,0], [ylim[0]+0.05]*24, "bv") #projection class2 to x-axis
    pl.plot([xlim[0]+0.03]*24, Dxn[L==-1,1], "b<")  #projection class2 to y-axis
    
    pl.plot(margin_xs, margin_ys, "g--", label="SVM Linear Separator")
    pl.xlim(xlim)
    pl.ylim(ylim)
    pl.xlabel(gl[0])
    pl.ylabel(gl[1])
    pl.legend()
    pl.title('Bivariate Separation of Samples using HLA-E and IL17RC')  
    pl.show()
    
def table_8_supplemental_all_pairs_classify():
    '''
    generates the classification rates for all pairs
    of the identified influenza pathways
    '''    
    resM = sp.zeros((12,12))
    
    #fill in the diagonal
    rc1,_,_ = all_pathway_classification()
    for i in range(12):
        resM[i,i] = round(rc1[i],3)
        
    #fill in the upper triangle
    rc2 = all_pathway_pairs_classification()
    for i in range(1,13):
        for j in range((i+1),13):
            #print i,j
            acc, _, _, _ = rc2[(i,j)]
            resM[i-1,j-1] = round(acc,3)
            #resM[j-1,i-1] = round(acc,3) #symmetric
        
    strM = (resM*100).astype("|S4")
    strM[ strM== "0.0" ] = "--"
        
    #collate row labels
    pw_dict = flu_genelists.PATHWAYS
    rowLabels = []    
    for i in range(12):
        r_str = "%d:~%s"%( (i+1), pw_dict[i+1][0])
        rowLabels.append(r_str)
            
    #build column-oriented structure for latex table output
    data_dict = {}
    data_dict['Pathway'] = ['1','2','3','4','5','6','7','8','9','10','11','12']
    for i in range(1,13):
        data_dict['%d'%i] = list( strM[:,(i-1)])
    
    cap = r"{\bf Classification accuracies of all pathway pairs.} "
    cap += "Values in bold face are those combinations that yield classifiers "
    cap += "with accuracies superior to the best accuracy of any single pathway. "
    cap += "Pathways are: %s"%(", ".join(rowLabels))
    try:
        asciitable.write(data_dict, sys.stdout, Writer = asciitable.Latex,
             col_align="|r|c|c|c|c|c|c|c|c|c|c|c|c|", latexdict = {'tabletype': 'table*', 
             'preamble': [r'\centering','\caption{%s}'%cap, r'\label{tbl:all_pathway_pairs}'],
             'data_start':r"\hline",
             'data_end':r"\hline", 'header_start':r"\hline" },
             names=['Pathway','1','2','3','4','5','6','7','8','9','10','11','12'])
    except:
        print "Error calling asciitable module for output of LaTeX."
  
    return resM

def figure_11_generanks():
    '''
    Generates supplemental figure showing
    the distribution of gene ranks for the top pathway
    and pathway pair classifiers.
    '''
    pl.subplot(1,2,1)
    _tmp = _pathways_generank()
    pl.subplot(1,2,2)
    _tmp = _pathway_pairs_generank()
    pl.show()


##############################################
#  Additional supporting functions
##############################################

def _pathways_generank():
    '''
    Helper function for figure 11
    '''
    (D,L,_,_) = ifr.load_flu_mat()
    G = ifr.load_gene_ids(short_name_only=True)
    dat = (D,L,_,G)
    pw_ranks = []
    data_vectors = []
    for k in ifr.flu_genelists.PATHWAYS.keys():
        pw, gl = ifr.flu_genelists.PATHWAYS[k]
        print pw
        rc = getGeneRanks(gl, dat)
        rnks = zip(*rc)[0]
        data_vectors.append(sp.array(rnks))
        print "Mean rank: %3.2f"%(sp.mean(sp.array(rnks)))
        pw_ranks.append( (pw,rnks) )
        print rnks
        
    pl.boxplot(data_vectors)
    pl.ylabel("Univariate Gene Rank (T-Test)")
    pl.xlabel("Discriminative Pathway")
    pl.ylim((0,12000))
    return pw_ranks
    
def _pathway_pairs_generank():
    '''
    Helper function for figure 11
    '''
    (D,L,_,_) = ifr.load_flu_mat()
    G = ifr.load_gene_ids(short_name_only=True)
    dat = (D,L,_,G)
    pwp_ranks = []
    data_vectors = []
    rc = table_3_pathway_pairs()
    pw_pairs = rc['Pathway Pair']
    d2 = dict( ifr.flu_genelists.PATHWAYS.values())
    for pwp in pw_pairs:
        p1,p2 = pwp.split(' + ')
        gl1 = d2[p1]
        gl2 = d2[p2]
        gl = gl1+gl2
        
        rc = getGeneRanks(gl,dat)
        rnks = zip(*rc)[0]
        data_vectors.append( sp.array(rnks) )
        print "Mean rank: %3.2f"%(sp.mean(sp.array(rnks)))
        pwp_ranks.append( (pwp,rnks) )
        print rnks

    pl.boxplot(data_vectors)
    pl.ylabel("Univariate Gene Rank (T-Test)")
    pl.xlabel("Discriminative Pathway Pair")
    pl.ylim((0,12000))
    return pwp_ranks

def gene_expressions_BEnet50():
    '''
    Generates the box plot of the gene expression levels for Chen et al 2011, Bayesian Elastic Net
    gene selection.
    '''
    (D,L,_,_) = ifr.load_flu_mat()
    G = ifr.load_gene_ids(short_name_only=True)
    gl = flu_genelists.BENET50_GENELIST
    idxs = [G.index(x) for x in gl]
    ifr.gene_expression_boxplots(D, L, gl , idxs, zscore=True, titlestr="Bayesian ENet Top 50")
    pl.show()
    

def gene_expressions_SSVM(genelist=IFR_INFLUENZA_FILE_GENES, itr=1):
    '''
    Generates the box plot of the gene expression levels for SSVM Iterative Feature Removal
    for a specified iteration (use one-based indexing to specify iteration).
    '''
    if itr < 1:
        raise ValueError("Iteration must be specified as >= 1.")
    (D,L,_,_) = ifr.load_flu_mat()
    G = ifr.load_gene_ids(short_name_only=True)
    removals = ifr.ifr_load_removals_from_file(genelist)  #list of lists, one per iteration
    gl = removals[itr-1]
    idxs = [G.index(x) for x in gl]
    ifr.gene_expression_boxplots(D, L, gl , idxs, zscore=True, titlestr="IFR Iteration %d"%itr)
    pl.show()

def gene_expressions_pathway(genelist, pathway_name,class_labels=('class1','class2'),
                             zscore=True, show=True, dat=None):
    '''
    Generates box plots of gene expressions for an arbitrary list of genes, with the intention of
    being able to list genes discovered from a known pathway, or t-test ranking, or any other
    selection criteria.
    @param genelist: A list of gene names
    @param pathway_name: A string describing what this set is
    @param class_labels: Tuple ('class1', 'class2'), where you should replace the
    strings with the appropriate class names for your data
    @param zscore: If true, the expression levels will be mean-centered and scaled to
    unit deviation. If false, the 'raw' expression levels will be used, but this can
    cause problems if the scales from one gene to the next vary widely.
    @param show: If true, the figure will be immediately drawn. If false, you can
    further annotate the figure with other matplotlib functions before display.
    @param dat: A tuple (D,L), where D is the data matrix having samples in rows
    and L is the label vector. If None, the H3N2 flu data will be used.
    '''
    if dat is None:
        (D,L,_,_) = ifr.load_flu_mat()
    else:
        (D,L) = dat
        
    G = ifr.load_gene_ids(short_name_only=True)
    idxs = [G.index(x) for x in genelist]
    ifr.gene_expression_boxplots(D, L, genelist , idxs, zscore=zscore, titlestr=pathway_name,
                                 class_labels=class_labels)
    if show: pl.show()

def getGeneRanks(gl, dat):
    '''
    For a given genelist, gl (names), and a data set
    tuple (D,L,_,G), returns the ranking according to t-test
    of each gene and the t statistic score.
    @param gl: The gene list, as gene names
    @param dat: A tuple, such as that returned by load_XYZ_dat()
    functions, (D,L,_Id,G), where the _Id field is ignored.
    '''
    G = dat[-1]
    gx = [ G.index(g) for g in gl ]
    return getGeneRanks2(gx, dat)
  
def getGeneRanks2(gx, dat):
    '''
    As per get GeneRanks, but uses as input a list of feature indices
    instead of gene names. 
    '''
    (D,L,_,G) = dat
    y = ifr.ttest_features(D, L, return_mat=False)
    res = [ (i,G[idx],tscore) for i,(tscore,_,idx) in enumerate(y) if idx in gx]
    return res
    
    
def genelist_iterations(gl, removals=None, fn=IFR_INFLUENZA_FILE_GENES):
    '''
    Output the iteration number for each gene in gl
    @param gl: The gene list to analyze
    @param removals: If None, the iterative removals will be loaded
    from the file specified by the fn parameter. Otherwise, provide
    the list-of-lists structure with the removal data.
    @param fn: If removals is None, then this is the file to load
    the data from.
    '''
    if removals is None:
        removals = ifr.ifr_load_removals_from_file(fn)  #list of lists, one per iteration
    itlist = []
    
    for g in gl:
        for i,rmvList in enumerate(removals):
            if g in rmvList: itlist.append(i)
            continue
    
    return itlist

def _load_stored_IFR(which):
    '''
    Function which can create an IFR object from
    csv files that store the genes removed per iteration.
    This is helpful, for example, when the SSVM version of
    IFR was run using matlab and the data saved to csv files.
    This function can load the csv files, instantiate an
    IFR object, which is then used in analysis to create
    figures, etc.
    @param which: Specify either 'flu' or 'lung', as these
    are the only two csv data files available to load.
    '''
    assert which in ('flu','lung')
    
    if which == 'flu':
        (dat1, dat2) = load_flu_data()
        genelist = IFR_INFLUENZA_FILE_GENES 
        acclist = IFR_INFLUENZA_FILE_ACC
    else:
        (dat1, dat2) = load_lung_data()
        genelist = IFR_LUNG_FILE_GENES 
        acclist = IFR_LUNG_FILE_ACC
           
    (D,L,_,_) = dat1
    (D2,L2,_,_) = dat2
    
    removals = ifr.ifr_load_removals_from_file(genelist)  #list of lists, one per iteration    
    A = sp.loadtxt(acclist)  #accuracy data as a matrix
    
    ifrx = ifr.IFR( (D,L), (D2,L2))
    ifrx.removed = removals
    ifrx.train_acc_list = list(A[:,3])
    ifrx.test_acc_list = list(A[:,4])
    
    return ifrx


# New functions to support additional analysis
def load_flu_data_11_14():
    '''
    Loads the H3N2 data and the H1N1 using the
    time intervals 11,12,13, and 14, which are
    better representative of "peak symptom" than
    the interval 14,15,16.
    H3N2 data can be used for training, validation,
    and model selection, with H1N1 witheld to test
    the best classifiers selected in validation.
    @return ( (D,L,Id,G), (D2,L2,Id2,G) ) where D,L,Id,G 
    are for training/validation and D2,L2,Id2,G are
    for testing.
    '''
    (D,L,_,Id) = ifr.load_flu_mat(times=[11,12,13,14])
    (D2,L2,_,Id2) = ifr.load_flu_mat(matfile='H1N1_dat.mat', times=[11,12,13,14], id_filter=[1,3,5,10,18])
    G = ifr.load_gene_ids(short_name_only=True)
    return ((D,L,Id,G), (D2,L2,Id2,G))

def load_flu_data():
    '''
    Loads the H3N2 data and the H1N1 using the
    time intervals 14,15,16, which are those used
    by Duke in their papers.
    H3N2 data can be used for training, validation,
    and model selection, with H1N1 withheld to test
    the best classifiers selected in validation.
    @return ( (D,L,Id,G), (D2,L2,Id2,G) ) where D,L,Id,G 
    are for training/validation and D2,L2,Id2,G are
    for testing.
    '''
    (D,L,_,Id) = ifr.load_flu_mat(times=[14,15,16])
    (D2,L2,_,Id2) = ifr.load_flu_mat(matfile='H1N1_dat.mat', times=[14,15,16], id_filter=[1,3,5,10,18])
    G = ifr.load_gene_ids(short_name_only=True)
    return ((D,L,Id,G), (D2,L2,Id2,G))

def load_lung_data():
    '''
    Loads the lung cancer data, and translates the
    feature names from affy ids to gene names using a stored
    dictionary that maps from one to the other.
    @return ( (D,L,Id,G), (D2,L2,Id2,G) ) where D,L,Id,G 
    are for training/validation and D2,L2,Id2,G are
    for testing.
    '''
    (X1,X2) = ifr.load_lung_cancer_dat()
    D = X1.data
    L = sp.array(X1.labelsAsIntegers())
    D2 = X2.data
    L2 = sp.array(X2.labelsAsIntegers())
    
    affy_dict = ifr.load_lung_cancer_dict()
    affy_ids = X1.fields
    G = [ affy_dict[x][2] for x in affy_ids]
    
    return ( (D,L,None,G), (D2,L2,None,G))

def load_prostate_data():
    '''
    Loads the prostate data, and converts the
    feature names to gene names, using the same
    affy_id to gene_id that was created for the
    lung data.
    @return ( (D,L,Id,G), (D2,L2,Id2,G) ) where D,L,Id,G 
    are for training/validation and D2,L2,Id2,G are
    for testing.
    '''
    (X1,X2) = ifr.load_prostate_dat()
    D = X1.data
    L = sp.array(X1.labelsAsIntegers())
    D2 = X2.data
    L2 = sp.array(X2.labelsAsIntegers())
    
    affy_dict = ifr.load_lung_cancer_dict()
    affy_ids = X1.fields
    G = [ affy_dict[x][2] for x in affy_ids]
    
    return ( (D,L,None,G), (D2,L2,None,G))

def load_bcell_data(stored_split=True):
    '''
    Loads bcell lymphoma data, which doesn't have a pre-defined
    train/test split. So we randomly pick 25% to be the withheld
    test partition. Maps the feature names to gene names using
    a stored dictionary file.
    @param stored_split: If True, then we have already stored
    a copy of this data with a particular train/test split, and to
    reproduce results, this data should be loaded. If False, a new
    split will be produced.
    @return ( (D,L,Id,G), (D2,L2,Id2,G) ) where D,L,Id,G 
    are for training/validation and D2,L2,Id2,G are
    for testing.
    '''
    if stored_split:
        fn = os.path.join(ifr.RESULTS_DIR, 'discriminative pathways', 'bcell', 'dat1_dat2.p')
        (dat1,dat2) = cPickle.load(open(fn,'rb'))
    else:
        X = ifr.load_bcell_lymphoma_dat()
        gdict = ifr.load_bcell_lymphoma_dict()
        G = [ gdict[f] for f in X.fields]
        
        D = X.data
        L = sp.array(X.labelsAsIntegers())
        
        (D1,L1,D2,L2,_) = ifr.random_partition(D, L, pctTrain=0.75, proportional_labels=True)
        dat1 = (D1,L1,None,G)
        dat2 = (D2,L2,None,G)
    
    return (dat1, dat2)

##############################################
# The following functions support the
# multiple trials (cross validation) results
# using IFR and automated pathway classification
# using GO annotations.
##############################################

def mt_bcell():
    '''
    Performs multiple-trial (50) cross-validation of pathway and pathway pair
    classifiers using the bcell lymphoma cancer data set. This is a top-level function that
    just sets some parameters and calls mt_script1 and mt_script2 and saves
    the results.
    '''
    base_dir = os.path.join(ifr.RESULTS_DIR,"discriminative pathways","bcell")
    fn1 = os.path.join(base_dir, "bcell_script1_res.p")
    (dat1, dat2) = load_bcell_data(stored_split=True)
    (pw_res, pwp_res, pc_list) = mt_script1(dat1, dat2, "bcell", t1=0.7, t2=0.75, depth=5,
                      min_pw_count=1, min_pw_size=5, common_scale=False,
                      numIter=20, save_as=fn1, save_pclist=True)
    
    d = mt_script2(pc_list, dat1, dat2, dataset='Bcell Lymphoma', common_scale=False)
    fn2 = os.path.join(base_dir, "bcell_script2_res.p")
    cPickle.dump(d, open(fn2,"wb"))
    return (pw_res, pwp_res, pc_list, d)
        

def mt_prostate():
    '''
    Performs multiple-trial (50) cross-validation of pathway and pathway pair
    classifiers using the prostate cancer data set. This is a top-level function that
    just sets some parameters and calls mt_script1 and mt_script2 and saves
    the results.
    '''
    base_dir = os.path.join(ifr.RESULTS_DIR,"discriminative pathways","prostate")
    fn1 = os.path.join(base_dir, "prostate_script1_res.p")
    (dat1, dat2) = load_prostate_data()
    (pw_res, pwp_res, pc_list) = mt_script1(dat1, dat2, "prostate", t1=0.7, t2=0.75, depth=5,
                      min_pw_count=1, min_pw_size=5, common_scale=False,
                      numIter=30, save_as=fn1, save_pclist=True)
    
    d = mt_script2(pc_list, dat1, dat2, dataset='Prostate Cancer', common_scale=False)
    fn2 = os.path.join(base_dir, "prostate_script2_res.p")
    cPickle.dump(d, open(fn2,"wb"))
    return (pw_res, pwp_res, pc_list, d)

def mt_lung():
    '''
    Performs multiple-trial (50) cross-validation of pathway and pathway pair
    classifiers using the lung cancer data set. This is a top-level function that
    just sets some parameters and calls mt_script1 and mt_script2 and saves
    the results.
    '''
    base_dir = os.path.join(ifr.RESULTS_DIR,"discriminative pathways","lung")
    fn1 = os.path.join(base_dir, "lung_script1_res.p")
    (dat1, dat2) = load_lung_data()
    (pw_res, pwp_res, pc_list) = mt_script1(dat1, dat2, "lung", t1=0.8, t2=0.95, depth=5,
                      min_pw_count=1, min_pw_size=5, common_scale=True,
                      numIter=30, save_as=fn1, save_pclist=True)
    
    d = mt_script2(pc_list, dat1, dat2, dataset='Lung Cancer', common_scale=True)
    fn2 = os.path.join(base_dir, "lung_script2_res.p")
    cPickle.dump(d, open(fn2,"wb"))
    return (pw_res, pwp_res, pc_list, d)

def mt_flu():
    '''
    Performs multiple-trial (50) cross-validation of pathway and pathway pair
    classifiers using the influenza 11-14 data set. This is a top-level function that
    just sets some parameters and calls mt_script1 and mt_script2 and saves
    the results.
    '''
    base_dir = os.path.join(ifr.RESULTS_DIR,"discriminative pathways","influenza")
    fn1 = os.path.join(base_dir, "flu_script1_res.p")
    (dat1, dat2) = load_flu_data_11_14()
    (pw_res, pwp_res, pc_list) = mt_script1(dat1, dat2, "flu", t1=0.8, t2=0.95, depth=5,
                      min_pw_count=1, min_pw_size=5, common_scale=False,
                      numIter=30, save_as=fn1, save_pclist=True)
    
    d = mt_script2(pc_list, dat1, dat2, dataset='Influenza 11-14', common_scale=False)
    fn2 = os.path.join(base_dir, "flu_script2_res.p")
    cPickle.dump(d, open(fn2,"wb"))
    return (pw_res, pwp_res, pc_list, d)
    
###########################################
# Core code to support the top-level
# functions for multiple trial (mt)
# cross validation of the pathway
# and pathway-pair classification
# results
###########################################
    
def mt_script1(dat1, dat2, ds_name, t1=0.9, t2=0.95, pc_list=None, depth=5, min_pw_count=1,
               min_pw_size=5, common_scale=False, numIter=30, save_as=None, save_pclist=True):
    '''
    Automates the process of multiple trial pathway classification analysis
    @param dat1: A tuple (D,L,Id,G) which defines the data used for cross-validated training
    (including feature selection) of the pathway classifiers. D is the data matrix, samples
    in rows, L is the integer label vector, Id can be None or an integer vector indicating
    the subject id of the samples (which allows the random partitioning to prevent subject
    overlap), and G is a list of the gene names, in the column (feature) order of D.
    @param dat2: As above, but this data is the withheld test data which is used after
    the top pathways were selected via x-validation.
    @param ds_name: A string to represent the data set in the result output, like "Influenza"
    or "Lung Cancer"
    @param t1: The accuracy threshold for picking top pathway classifiers
    @param t2: The accuracy threshold for picking top pathway pairs
    @param pc_list: If you have a pre-computed list of pathway classifier objects, such as
    output by multiple_trials_pathway_classification, you can skip that aspect of this
    script.
    @param min_pw_count: The minimum number of trials that a gene must appear in order to be
    considered part of the pathway. The gene must show up in at least this many trials to be included,
    the default is 1, meaning that even if the gene only appeared once in a single trial, it will
    still be included in the discriminative pathway.
    @param min_pw_size: Pathways must have at least this number of genes in order to be considered.
    @param save_as: Optional. Provide a string representing a pickle file name/path where
    the results should be stored.
    @param save_pclist: Ignored unless save_as is not None. If true, the pc_list data structure
    is saved as part of the results, if false the pc_list data is not saved. The pc_list data takes
    up a large amount of disk space (perhaps 300-500 megabytes), and, especially if you pre-computed
    the pc_list and passed it into this function, you don't want to store it again.
    @return: (pw_res, pwp_res, pc_list)
    '''
    (D,L,_,G) = dat1
    (D2,L2,_,_) = dat2
    
    if pc_list is None:
        print "Step 1: Multiple trial IFR and Pathway Classification"
        pc_list,_ = multiple_trials_pathway_classification(dat1, nTrials=50, depth=depth, common_scale=True,
                                                           numIter=numIter)
        #Note: the top-level common_scale parameter is not used here, it is used for normalizing the test
        # data...whether the test data should be normalized using the mean/std of training, or be 'self-centered'.
        # The latter is often better when the test data is from a different collection than the training/validation data.
     
     
    print "Step 2: Finding the best average pathway classifiers"
    A, pw = mt_get_A_matrix(pc_list)
    #(Am,As,_Ns) = ifr.gappy_means_and_stdvs(A, sentinel=0.0) #<-gappy version
    (Am,As,_Ns) = ifr.gappy_means_and_stdvs(A, sentinel=-1.0) #<-sentinel doesn't exist, so non-gappy version
    
    best_pws = [ (p,a,s) for (p,a,s) in zip(pw,Am,As) if a >= t1 ]    
    print "There are %d pathways with %3.2f or better mean accuracy"%(len(best_pws), t1)
    
    print "Step 3: Evaluating best pathways on withheld test data"
    pw_res = []
    for (p,a,s) in best_pws:
        rc = mt_pathway_profile(pc_list, p, show=False)
        gl = [ g for (c,g) in rc[0] if c >= min_pw_count ]  #gene list
        if len(gl) < min_pw_size:
            print "Pathway %s only has %d genes, less than min_pw_size. It will be skipped."%(p,len(gl))
            continue
        
        res = ifr.pathway_classification(gl, feature_names=G, traindat=(D,L), testdat=(D2,L2),
                                         common_scale=common_scale, verbose=False)
        test_acc = res[0]
        ci = 1.96 * s / math.sqrt( len(pc_list) )
        print "%s \t x-val:%3.1f +- %3.1f \t test: %3.1f"%(p, 100*a, 100*ci, 100*test_acc)
        pw_res.append( (p, a, s, test_acc) )
        
    print "Step 4: Finding best pathway pair classifiers"
    AP, pwp = mt_get_AP_matrix(pc_list, thresh=t1) #when compiling pairs, use only pathways from prev step
    APm = AP.mean(axis=0)
    APs = AP.std(axis=0)
    best_pps = [ (pp,a,s) for (pp,a,s) in zip(pwp,APm,APs) if a >= t2 ] #higher cutoff as these tend to be more accurate
    print "There are %d pathway pairs with %3.2f or better mean accuracy"%(len(best_pps), t2)
    
    print "Step 5: Evaluating best pathway pairs on withheld test data"
    pwp_res = []
    for ((p1,p2),a,s) in best_pps:
        pp_str = "%s + %s"%(p1[:18],p2[:18])
        rc1 = mt_pathway_profile(pc_list, p1, show=False)
        rc2 = mt_pathway_profile(pc_list, p2, show=False)
        
        gl1 = [ g for (c,g) in rc1[0] if c >= min_pw_count ]  #gene list
        gl2 = [ g for (c,g) in rc2[0] if c >= min_pw_count ]  #gene list
        gl3 = list( set(gl1).union(set(gl2)) )
        
        if len(gl1) < min_pw_size or len(gl2) < min_pw_size: continue
                
        if len(gl3) < min_pw_size:
            print "Pathway pair %s only has %d genes, less than min_pw_size. It will be skipped."%(pp_str,len(gl3))
            continue
                
        res = ifr.pathway_classification(gl3, feature_names=G, traindat=(D,L), testdat=(D2,L2),
                                         common_scale=common_scale, verbose=False)
        test_acc = res[0]
        ci = 1.96 * s / math.sqrt( len(pc_list) )
        print "%s \t x-val:%3.1f +- %3.1f \t test: %3.1f"%(pp_str, 100*a, 100*ci, 100*test_acc)
        pwp_res.append( ((p1,p2), a, s, test_acc) )    
    
    res =  (pw_res, pwp_res, pc_list)
    if not save_as is None:
        if save_pclist:
            cPickle.dump(res, open(save_as,"wb"))
        else:
            cPickle.dump((pw_res,pwp_res,None), open(save_as,"wb"))
        
    return res

def mt_script2( pc_list, dat1, dat2, dataset='influenza', common_scale=False):
    '''
    This script summarizes the test results and gene ranking of features
    involved with the best svm classifier on the data set (sparse svm), compared to
    the best pathway and pathway pair classifier.
    
    The results are intended to show that deeper feature mining yields classifiers
    with accuracy on par of the 'optimal' classifier, but using sets of genes with
    much lower t-test gene ranking (i.e., less obvious features).
    '''
    d = {}
    print "Cross-validation testing using sparse svm classifier..."
    (_,acc_list,_,gr_svm) = cross_validation_score(dat1, dat2, r=50, common_scale=common_scale)
    best_svm_test_acc = acc_list[0]
    
    d['dataset'] = dataset
    d['generank_svm'] = gr_svm
    d['test_acc_svm'] = best_svm_test_acc
    
    print "Computing best pathway and pathway pair classifiers..."
    rc = best_pw_pwp_generanks(pc_list, dat1, dat2, common_scale=common_scale, min_pw_count=2,
                               min_pw_size=5)
    ( (pw, pw_test_acc, _, gr_pw), ((p1,p2), pwp_test_acc, _, gr_pp)) = rc
    
    d['pathway'] = pw
    d['pair'] = "%s + %s"%(p1,p2)
    d['generank_pathway'] = gr_pw
    d['generank_pair'] = gr_pp
    d['test_acc_pathway'] = pw_test_acc
    d['test_acc_pair'] = pwp_test_acc
    
    return d

def multiple_trials_pathway_classification(dat, nTrials=10, depth=5, numIter=30, common_scale=True):
    '''
    Performs the entire process of IFR feature selection +
    ontology analysis for pathways + computing pathway and pathway pair
    classification accuracies. For each trial, a random 75/25% subject-based
    splitting is used for the data. Subject splits are proportional between
    the classes, so as to have approximately balanced classes in training and
    validation data.
    @param dat: Tuple (D,L,Id,G) where D is the data matrix, L is the label vector,
    and Id is the list of subject ids. If Id is None, then the data is assumed to
    have no subject overlap. Otherwise, the ids are used to perform subject-based
    random partitioning of the data to ensure no subject overlap. G is the list of
    gene names in order of the columns of D.
    @param nTrials: How many pathway classifiers should be constructed using
    random partitioning of D/L.
    @param depth: The ontology depth used for annotations of the pathways.
    @return: (pc_list, id_list), where pc_list is a list of pathway classifier objects,
    one per trial, and id_list is the associated set of validation subject
    identifiers used per iteration.
    '''
    (D,L,Id,G) = dat
    
    pc_list = []
    id_list = []
    
    for i in range(nTrials):
        print "========================"
        print "=       TRIAL %s       ="%str(i+1).zfill(2)
        print "========================"
        if Id is None:
            (trainD, trainL, validD, validL, _) = ifr.random_partition(D, L, pctTrain=0.75,
                                                                            proportional_labels=True)
        else:
            (trainD, trainL, validD, validL, ids) = ifr.subject_partition_proportional(D, L, Id, 0.25)
            id_list.append(ids)
        pc = ifr.Pathway_Classifier( (trainD,trainL), (validD,validL), G, numIter=numIter)
        pc.initAnalysis( IFR_kwargs={'common_scale':common_scale, 'C':0.5},
                         IFRA_kwargs={'depth':depth} )
        pc.computePathwayClassifiers(min_size=1)
        pc.computePathwayPairsClassifiers()
        pc_list.append(pc)
        
    return (pc_list, id_list)
    
def mt_get_A_matrix(pc_list):
    '''
    From the results of multiple_trials (mt) pathway classification,
    compute the matrix A where the entry A[i,j] is
    the pathway classification accuracy for trial i, pathway j.
    @param pc_list: The list of pathway classifiers computed via 
    multiple_trials_pathway_classification() function.
    @return: A, all_pathways
    '''
    tmp = [ pc.pathways.keys() for pc in pc_list ]
    all_pathways = sorted( set([item for sublist in tmp for item in sublist]))
    p = len(all_pathways)
    N = len(pc_list)
    A = sp.zeros((N,p))
    
    for row,pc in enumerate(pc_list):
        for (acc, pw) in pc.getTopPathways(1000):
            col = all_pathways.index(pw)
            A[row,col] = acc
      
    return A, all_pathways
    
def mt_get_AP_matrix(pc_list, thresh=0.6, verbose=True):
    '''
    From the results of multiple_trials (mt) pathway classification,
    compute the matrix AP where the entry AP[i,j] is
    the pathway PAIR classification accuracy for trial i, pathway PAIR j.
    @param pc_list: The list of pathway classifiers computed via 
    multiple_trials_pathway_classification() function.
    @return: AP, all_pathway_pairs
    '''
    A, pw = mt_get_A_matrix(pc_list)
    
    #sub-select pw to those with threshold accuracy
    pw_idxs = sp.nonzero( A.mean(axis=0) >= thresh )[0]
    pw2 = sorted( [ pw[i] for i in pw_idxs ] )
    
    #Method 1
    #get list of genes for each pathway in pw2
    #perform all_pairs classification with aggregate lists
    
    #Method 2
    #from each pc in pc_list, find pathway_pair score
    # stack in a matrix with zero values where appropriate
    N = len(pw2)
    m = ( N*(N-1) ) /2
    AP = sp.zeros((len(pc_list),m))
    pw_pairs_gen = it.combinations(pw2, 2) # all combinations 2-at-a-time, generator
    pwps = []    
    for col,(p1,p2) in enumerate( pw_pairs_gen ):
        pwps.append( (p1,p2) )
        key = "%s + %s"%(p1,p2)
        key_alt = "%s + %s"%(p2,p1)
        if verbose: print "Compiling pathway pair: %s"%key
        for row,pc in enumerate(pc_list):
            d = {pwp:acc for (acc,pwp,_) in pc.pathway_pair_accuracies.values() }
            if key in d:
                AP[row,col] = d[key]                
            elif key_alt in d:
                AP[row,col] = d[key_alt]     
            
    return AP, pwps
    
    
    
def mt_pathway_profile(pc_list, pw, show=False):
    '''
    Shows the profile of a given pathway. The profile
    is composed of the genes in the pathway (x-axis)
    by the number of times over all the trials that
    gene was selected by IFR (y-axis).
    @param pc_list: A list of pathway classifier result
    objects, one per trial, as output by multiple_trials_pathway_classification()
    @param pw: The pathway name to analyze
    @param show: If true, the profile will be presented as
    a figure. If false, just profile information will be
    returned.
    @return: Tuple: (gene_counts, pw_genes, percents)
    '''
    pw_genes = []
    for pc in pc_list:
        if pw in pc.pathways:
            pw_genes.append( pc.pathways[pw])
        else:
            pw_genes.append( [] )
            
    #all genes relating to this pathway over all trials
    all_pw_genes = sorted(set([ item for sublist in pw_genes for item in sublist]))
    N = float( len(all_pw_genes) )
    
    #for each trial, count the pct of all_pw_genes found in that trial
    pcts = [ len(gl)/N for gl in pw_genes ]
    
    #for each gene in all_pw_genes, find out how many trials it appears            
    gene_counts = [ (sum([(g in gl) for gl in pw_genes]),g) for g in all_pw_genes]
    
    if show: mt_plot_profile(pw, gene_counts, numtrials=len(pc_list))
    
    return (gene_counts, pw_genes, pcts)

def mt_plot_profile(pw, gene_counts,numtrials=20):
    '''
    Helper function to mt_pathway_profile that plots the profile
    generated.
    '''
    counts, genes = zip(*gene_counts)
    N = len(counts)
    pl.bar( range(N), counts )
    pl.xlim( (0,N+1) )
    pl.ylim( (0,numtrials))
    pl.ylabel('Number of trials (of %d) where pathway contains gene'%numtrials)
    pl.title('Pathway profile %s'%pw)
    pl.xticks(pl.arange(N)+0.5, genes, rotation='vertical', size='x-small',
              horizontalalignment='center')
    pl.subplots_adjust(top=0.95, bottom=0.12, left=0.05, right=0.95)
    pl.show()
    
def mt_plot_results(A, all_pathways, thresh=0.6, min_pw_size=1, title=None, show=True):
    '''
    Plots a bar chart with error bars for the result matrix A,
    from mt_get_A_matrix()
    @param A: The accuracy matrix from mt_get_A_matrix() function
    @param all_pathways: The list of all pathways used in A, also returned
    by mt_get_A_matrix()
    @param thresh: Set this value as the minimum mean accuracy a pathway
    must have to be shown in the bar chart. None = show all.
    @param show: If true, the plot will be displayed, else you have to
    call pl.show() at some point to see it. 
    '''
    Am = A.mean(axis=0)
    Amed = sp.median(A, axis=0)
    As = A.std(axis=0)
    Ae = ( As / math.sqrt(len(As)))
    apw = sp.array(all_pathways)
    
    if not thresh is None:
        mask = (Am >= thresh)
        X = Am[mask]
        E = Ae[mask]
        P = apw[mask]
        Md = Amed[mask]
        #A2 = A[:,mask]
    else:
        X = Am
        E = Ae
        P = apw
        Md = Amed
        #A2 = A
    
    M = len(X)
    pl.figure(1, figsize=(10,6) )
    pl.bar( pl.arange(M), X, yerr=E, color='#BFBFBF', ecolor='red', align='center') #, width=0.75)
    pl.plot(pl.arange(M), Md, 'rd')
    pl.xticks([])   
    pl.ylim([0,1.05]) 
    for ix,s in enumerate(P):
        pl.text( ix, 0.05, s, rotation='vertical', horizontalalignment='center',
                 verticalalignment='bottom', color='blue')
        pl.text( ix, X[ix], "%2.1f"%(100*X[ix]), rotation='horizontal', horizontalalignment='center',
                 verticalalignment='top', color='blue')
    pl.ylabel('Accuracy')
    pl.xlabel('Pathway')
    if title: pl.title(title)
    pl.subplots_adjust(right=0.95, left=0.10, bottom=0.05, top=0.95)
    #pl.figure(2)    
    #pl.ylim([0,1.2])
    #pl.boxplot(A2, vert=True, widths=0.5, hold=True )
        
    if show: pl.show()
    
def mt_plot_results2(A, all_pathways, thresh=0.6, show=True):
    '''
    Plots a bar chart with error bars. This variant computes the mean/std
    using only the non-zero values in A. Sometimes a pathway shows up in
    the Go analysis only for a few iterations, but when it does, it is highly
    discriminative. This allows those pathways to be visualized in the results,
    but makes no penalty for the fact that the pathway is often not found via IFRA.
    @param A: The accuracy matrix with missing values as 0.0
    @param all_pathways: The list of all pathways in same order as X
    @param thresh: Set this value as the minimum mean accuracy a pathway
    must have to be shown in the bar chart. None = show all.
    @param show: If true, the plot will be displayed, else you have to
    call pl.show() at some point to see it. 
    '''
    
    #first we have to compute the mean and std for each pathway,
    # where we ignore the times when the pathway wasn't present in
    # an ifr result.
    (Mu,Sig,Ns) = ifr.gappy_means_and_stdvs(A, sentinel=0.0)
    
    tmp = [ (m, s/math.sqrt(n)) for (m,s,n) in zip(Mu,Sig,Ns) ] # (mean, std error) per pathway
            
    T = sp.array(tmp)
    X = T[:,0]
    E = T[:,1]    
    apw = sp.array(all_pathways)
    
    if not thresh is None:
        mask = (X >= thresh)
        X = X[mask]
        E = E[mask]
        P = apw[X >=thresh]
    else:
        P = apw
    
    M = len(X)
    pl.bar( pl.arange(M), X, yerr=E, color='#BFBFBF', ecolor='red')
    pl.xticks([])
    for ix,s in enumerate(P):
        pl.text( ix+0.4, 0.05, s, rotation='vertical', horizontalalignment='center',
                 verticalalignment='bottom', color='blue')
        pl.text( ix+0.4, X[ix], "%2.1f"%(100*X[ix]), rotation='horizontal', horizontalalignment='center',
                 verticalalignment='top', color='blue')
        
    if show: pl.show()
    
def mt_best_validation_score( pw_res, K=1 ):
    '''
    Returns the K-best pathways/pathway_pairs with their validation scores.
    @param pw_res: The pathway/pair results object, as returned by
    mt_script1(...)
    @param K: The top K will be returned.
    '''
    tmp = [ (va,s,ta,p) for (p,va,s,ta) in pw_res]
    return sorted(tmp, reverse=True)[0:K]

def mt_topK_average(pw_res, K=5):
    '''
    Averages the validation scores for the top K pathway or pathway pair classifiers.
    @return: The average classification accuracy of the top K pathway/pathway-pairs in the
    results object.
    '''
    rc = mt_best_validation_score(pw_res, K=K)
    scores = [s for (_,_,s,_) in rc]
    return sp.mean(scores)


    

def cross_validation_score(dat1, dat2, r=50, common_scale=True):
    '''
    Perform r cross-validation trials using randomized partitioning on dat1,
    then choose the best model to predict on dat2
    '''
    (D,L,_,_) = dat1
    L = sp.array(L)
    Dn,mu,sd = ifr.normalize_data(D)
    (D2,L2,_,_) = dat2
    
    f = ifr.random_partition_idxs
    cvg = ifr.CV_Generator(r, f, L, frac=0.25, proportional_labels=True)
    
    res_list = []
    
    for (t1, t2) in cvg:
        #for each trial with train/validation split t1/t2
        Dtrn = Dn[t1,:]
        Ltrn = L[t1]
        Dval = Dn[t2,:]
        Lval = L[t2]
        
        (acc,_,clf,_,_) = ifr.svm_engine((Dtrn,Ltrn), (Dval,Lval), common_scale=True,
                             no_normalization=True, loss='L2', penalty='L1', C=1.0 )   #L1 penalty induces sparseness
        
        res_list.append( (acc,clf) )
        
    res_list = sorted(res_list, reverse=True)  #best accuracy first
    A = sp.array([ a for (a,_) in res_list])
    Am = A.mean()
    As = A.std()
    print "Mean/Std of %d trials: %3.2f / %3.2f"%(r, Am*100, As*100)
    
    (best_a, best_clf) = res_list[0]  #best classifier and associate accuracy from validation trials
    print "Best classifier from validation trails scored %3.2f."%(best_a*100)
    
    print "Applying top 5 classifiers to withheld test data."
    if common_scale:
        D2n,_,_ = ifr.normalize_data(D2, mu, sd)  #normalize and center to training distribution
    else:
        D2n,_,_ = ifr.normalize_data(D2) #normalize/scale to own testing distribution
        
    res2 = []
    for i in range(5):
        (_,clf) = res_list[i]
        preds = clf.predict(D2n)
        score = float( sum(preds==L2))/ len(L2)
        
        res2.append(score)
    
    print "Best model test score: %3.2f"%(100*res2[0])
    print "Top 5 Avg test score: %3.2f"%(100*sp.mean(res2))

    print "T-Test ranking of genes in best model:"
    factors = sp.nonzero(best_clf.coef_[0,:])[0]
    gr_info = getGeneRanks2(factors, dat1)
    gr = [ x for (x,_,_) in gr_info]
    print gr
        
    return res_list, res2, gr_info, gr

def best_pw_pwp_generanks(pc_list, dat1, dat2, min_pw_count=2, min_pw_size=5, common_scale=False):
    '''
    Get the gene ranks for the genes in the best pathway and
    pathway pair extracted from a pc_list of cross-validated
    pc classification results
    '''
    A, pw = mt_get_A_matrix(pc_list)
    Am = A.mean(axis=0)
    
    ranked_pws = sorted(zip(Am,pw), reverse=True)  
    #look through the ranked pathways, best first, and find
    # the first that meets the criteria of min number of genes
    success = False
    for (_,pw) in ranked_pws:
        rc = mt_pathway_profile(pc_list, pw, show=False)
        gl = [ g for (c,g) in rc[0] if c >= min_pw_count ]  #gene list
        if len(gl) >= min_pw_size:
            success = True
            break
        
    if not success:
        print "Error: No pathways meet selection criteria!"
        return None
    
    print "Pathway: \'%s\' is the best ranked pathway meeting selection criteria"%pw
    (D,L,_,G) = dat1
    (D2,L2,_,_) = dat2
    svm_res = ifr.pathway_classification(gl, feature_names=G, traindat=(D,L), testdat=(D2,L2),
                                    common_scale=common_scale)
    
    print "Test accuracy is %3.2f"%svm_res[0]
    factors = svm_res[1]
    gr_info = getGeneRanks2(factors, dat1)
    gr = [x for (x,_,_) in gr_info]
    
    #ditto for the best pathway pair
    print "Compiling pairs..."
    AP, pwp = mt_get_AP_matrix(pc_list, thresh=0.7)
    APm = AP.mean(axis=0)
    ranked_pairs = sorted(zip(APm,pwp))
    success = False
    for (_,(p1,p2)) in ranked_pairs:
        rc1 = mt_pathway_profile(pc_list, p1, show=False)
        gl1 = [ g for (c,g) in rc1[0] if c >= min_pw_count ]  #gene list
        rc2 = mt_pathway_profile(pc_list, p2, show=False)
        gl2 = [ g for (c,g) in rc2[0] if c >= min_pw_count ]  #gene list
        gl3 = list( set(gl1).union(set(gl2)) )
        
        if len(gl1) < min_pw_size or len(gl2) < min_pw_size: continue
        
        if len(gl3) >= min_pw_size:
            success = True
            break
    
    if not success:
        print "Error: No pathway pairs meet selection criteria!"
        return None
        
    pp = "%s + %s"%(p1,p2)
    print "Pathway Pair: \'%s\' is the best ranked pathway meeting selection criteria"%pp
    (D,L,_,G) = dat1
    (D2,L2,_,_) = dat2
    svm_res2 = ifr.pathway_classification(gl3, feature_names=G, traindat=(D,L), testdat=(D2,L2),
                                    common_scale=common_scale)
    print "Test accuracy is %3.2f"%svm_res2[0]
    factors2 = svm_res2[1]
    gr_info2 = getGeneRanks2(factors2, dat1)
    gr2 = [x for (x,_,_) in gr_info2]
    
    return ( (pw, svm_res[0], factors, gr), ((p1,p2), svm_res2[0], factors2, gr2) )

def mt_plot_generank_res( gr_res_list, show=True):
    num_plots = len(gr_res_list)
    
    count = 1
    for ix in range( num_plots):
        pl.subplot(1,num_plots,ix+1)
        if count == 1:
            pl.ylabel("T-Test Gene Rank")
            count += 1
        d = gr_res_list[ix]
        gr1 = d['generank_svm']
        gr2 = d['generank_pathway']
        gr3 = d['generank_pair']
        #maxy = max( gr1+gr2+gr3 )
        ds = d['dataset']
        print ds
        pl.boxplot( (gr1,gr2,gr3))
        pl.xticks( [1,2,3], ['SSVM','Pathway','Pair'], rotation=65)
        #pl.ylim([0,maxy*1.1])
        #pl.text(0.75,maxy*1.05,"%3.2f"%d['test_acc_svm'])
        #pl.text(1.75,maxy*1.05,"%3.2f"%d['test_acc_pathway'])
        #pl.text(2.75,maxy*1.05,"%3.2f"%d['test_acc_pair'])
        pl.title(ds)

    pl.subplots_adjust(left=0.08, right=0.95, wspace=0.3, bottom=0.15)
    if show: pl.show()

def ttest_table6_res():
    datasets = {'lung':'pc_results_lung_80_95_26Oct2013.p',
                'influenza':'pc_results_flu_80_95_26Oct2013.p',
                'prostate':'pc_results_prostate_70_75_26Oct2013.p',
                'bcell':'pc_results_bcell_70_75_iter20_26Oct2013.p'}
    res = {}
    for ds in datasets.keys():
        fn1 = os.path.join(ifr.RESULTS_DIR,'discriminative pathways',ds,datasets[ds])
        fn2 = os.path.join(ifr.RESULTS_DIR,'discriminative pathways',ds,'ssvm_xval_res.p')
        (pw_res, pwp_res, _) = cPickle.load(open(fn1,'rb'))
        d = cPickle.load(open(fn2,'rb'))
        #get list of top 5 test scores from each method
        s1 = [ s for (_,_,s,_) in mt_best_validation_score(pw_res,5)]
        s2 = [ s for (_,_,s,_) in mt_best_validation_score(pwp_res,5)]
        s3 = d['test_accuracy_top5']
        #perform two-tail welch ttest on all pairs
        x12 = ttest_ind(s1, s2, equal_var=False)
        x13 = ttest_ind(s1, s3, equal_var=False)
        x23 = ttest_ind(s2, s3, equal_var=False)
        res[ds] = (x12, x13, x23)
        
    return res

def export_mt_data_to_csv(dataset_str, pc_list, pw_res, pwp_res, out_dir, t1=.80):
    '''
    For multiple-trial data, we want to export the results
    of the best pathways (and pathway pairs) over all
    N trials to a text file.
    @param dataset_str: A string like "Lung Cancer"
    @param pc_list: Loaded from a saved results, as per mt_script1
    @param pw_res: Also can be found in the saved results from mt_script1,
    this is the validation scores for each pathway
    @param pwp_res: As above, but for pathway pairs
    @param out_dir: A directory where the output files will be written,
    we suggest a separate sub-dir for each dataset, as the file names
    will not be unique
    @param t1: The threshold on the accuracy of a pathway for it to
    be considered as part of a pathway pair...use the same as was used
    in mt_script1.
    @return: Returns nothing, but creates csv data files in the out_dir.
    There will be one for the 50-trials of the pathways and another for
    the pathway pairs.
    '''
    (_,_,_,p) = mt_best_validation_score(pw_res)[0]
    (_,_,_,pp)= mt_best_validation_score(pwp_res)[0]
    
    A, pw = mt_get_A_matrix(pc_list)
    pidx = pw.index(p)
    tmp = os.path.join(out_dir, 'pw_list.txt')
    with open(tmp,'w') as f:
        outtxt = "\n".join(pw)
        f.write(outtxt)
    
    AP, pwp = mt_get_AP_matrix(pc_list, thresh=t1, verbose=True)
    ppidx = pwp.index(pp)
    tmp = os.path.join(out_dir, 'pwp_list.txt')
    with open(tmp,'w') as f:
        outtxt = "\n".join([str(x) for x in pwp])
        f.write(outtxt)
    
    ofile = os.path.join(out_dir, 'pw_scores.csv' )
    header_str = "Pathway %d-Trial Cross Validation for %s data.\n"%(len(pc_list), dataset_str)
    header_str += "Best: %s (1-based index: %d)\n"%(p,pidx+1)
    header_str += "Remember to count rows after stripping header..."
    sp.savetxt(ofile, A.T, delimiter=",", fmt="%3.4f", header=header_str)
        
    ofile2 = os.path.join(out_dir, 'pwp_scores.csv')
    header_str2 = "Pathway Pair %d-Trial Cross Validation for %s data.\n"%(len(pc_list), dataset_str)
    header_str2 += "Best: %s (1-based index: %d)\n"%(str(pp), ppidx+1)
    header_str2 += "Remember to count rows after stripping header..."
    sp.savetxt(ofile2, AP.T, delimiter=",", fmt="%3.4f", header=header_str2)
    
def export_ssvm_xval_to_csv(data_str='lung'):
    base_dir = os.path.join(ifr.RESULTS_DIR, 'discriminative pathways', data_str)
    fn = 'ssvm_xval_res.p'
    d = cPickle.load( open(os.path.join(base_dir,fn), 'rb'))
    xval_res = d['validation_scores_50trial']
    header_str = "SSVM %d-Trial Cross Validation on %s data.\n"%(len(xval_res),data_str)  
    header_str += "Each row is the validation accuracy of 1 trial.\n"
    out_data = ["%3.4f"%x for x in xval_res]
    
    outfile = os.path.join(base_dir,'ssvm_xval.txt')
    with open(outfile,'w') as f:
        f.write(header_str)
        f.write("\n".join(out_data))
    
def two_sample_ttest(data_str='lung', pidx=158, ppidx=72):
    base_dir = os.path.join(ifr.RESULTS_DIR, 'discriminative pathways', data_str)
    fn = 'ssvm_xval_res.p'
    d = cPickle.load( open(os.path.join(base_dir,fn), 'rb'))
    C = sp.array(d['validation_scores_50trial'])
    
    fn1 = os.path.join(base_dir,'pw_scores.csv')
    M1 = sp.loadtxt(fn1, comments='#', delimiter=',')
    A = M1[pidx,:]
    
    fn2 = os.path.join(base_dir,'pwp_scores.csv')
    M2 = sp.loadtxt(fn2, comments='#', delimiter=',')
    B = M2[ppidx,:]
    
    print "Welch's Two sample t-test between PW and PWP"
    rc1 = ttest_ind(A,B, equal_var=False)
    print rc1
    
    print "Welch's Two sample t-test between PW and SSVM"
    rc2 = ttest_ind(A,C, equal_var=False)
    print rc2
    
    print "Welch's Two sample t-test between PWP and SSVM"
    rc3 = ttest_ind(B,C, equal_var=False)    
    print rc3
    
    return (rc1,rc2,rc3)
    
if __name__ == '__main__':
    pass
    
    
    