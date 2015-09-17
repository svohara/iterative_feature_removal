'''
Created on Sep 12, 2013
@author: Stephen O'Hara

Code to help with plotting distribution of gene expressions and other
datat related to this code base.
'''
import ifr
from ifr.common.utilities import ttest_features
import pylab as pl
import scipy as sp

def gene_expression_boxplots(D, L, genelist, idxlist, zscore=False,
                    titlestr='Distribution of Expression Levels for Selected Genes',
                    class_labels=['Class 1','Class 2']):
    '''
    Function to show the univariate separability of a set of genes. For each gene in the list,
    two box plots are drawn, one per class. The box plots show the distribution of expression
    levels of the gene for each class, which gives an indicator of the univariate separability.
    @param D: The data matrix, rows=samples, cols=genes/features
    @param L: Labels, list or vector indicating class membership. Two-class analysis only.
    @param genelist: The text names of the genes
    @param idxlist: The column indexes into D associated with each gene in genelist
    @param zscore: If true, the graph will use normalized values (mean-centered, scaled by std),
    also known as "z-scores".
    @param titlestr: The title for the resulting graph
    '''
    #check to make sure labels are valid...integers for 2 classes, -1 / +1
    assert len(set(L)) == 2
    binary_labels = list(set(L))
    #print binary_labels
    
    #get ttest scores
    _, TT = ttest_features(D, L, return_mat=True)
        
    if zscore:
        D,_,_ = ifr.normalize_data(D)
    
    N = len(idxlist)
    #the following list completion is a complex nesting of commands to handle a small number
    # of genes that are too long to conveniently display, such as AFFX-HUMISGF3A/M97935_3_at,
    # in which case "M97935_3_at" will be displayed. Simple names like "OAS1" should be unaffected.
    class1_data = D[:,idxlist][L==binary_labels[1],:]
    odds = range(1,2*N,2)  #class1 boxes will be at the odd ticks
    class2_data = D[:,idxlist][L==binary_labels[0],:]
    evens = range(2, (2*N)+1,2) #class2 boxes will be at the even ticks
    
    genelist2 = [ gg.split("/")[1][:12] if "/" in gg else gg[:12] for gg in genelist]
    pl.hold(True)
    
    #draw boxplots for class1 at the odd ticks
    rc= pl.boxplot( class1_data, sym="b+",patch_artist=True,
                    positions=odds)  #the values for class -1
    for p in rc['boxes']: p.set_color("b")
    
    #draw boxplots for class2 at the even ticks
    rc2=pl.boxplot( class2_data, sym="rx",patch_artist=True,
                    positions=evens)  #the values for class +1
    for p in rc2['boxes']: p.set_color("r")
        
    #draw light vertical lines to group the box plot pairs
    ax = pl.gca()
    yrange = ax.yaxis.get_majorticklocs()
    pl.vlines(sp.array([0]+evens)+0.5, yrange.min(), yrange.max(), color='lightgray', linestyles='solid',
              alpha=0.5)

    #draw the ttest-measure score for each feature at the top
    x_pos_list = sp.array(range(1,2*N+1,2))+0.5
    y_pos = 0.9*yrange.max()
    for i,idx in enumerate(idxlist):
        tscore = TT[idx,0]
        #pval = TT[idx,1]
        pl.text( x_pos_list[i], y_pos, "%2.1f"%tscore, color='blue', ha='center')
    
    #labels, titles, and tick marks
    pl.xlabel("Gene")
    pl.ylabel("Expression Level")
    pl.title(titlestr)
    pl.xticks( x_pos_list, genelist2, rotation='vertical')
    pl.xlim([0,2*N+1])

    #legend
    r1 = pl.Rectangle((0,0),1,1,fc='b')
    r2 = pl.Rectangle((0,0),1,1,fc='r')    
    pl.legend( [r1,r2],class_labels, loc=4)

    pl.subplots_adjust(bottom=0.2, left=0.05, right=0.95)
    pl.draw()