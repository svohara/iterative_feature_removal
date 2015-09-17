'''
Created on Jan 2, 2013
@author: Stephen O'Hara

Functions to assist in converting affymetrix probe_ids and other data elements
into standard ids that can be used to query public databases
'''
import os
import ifr
import cPickle

def get_geneids_from_affy( affy_id_list, affy_file=None ):
    '''
    Returns a dictionary mapping affy probe ids to
    the tuple (genebank,unigene,symbol) given an input list
    of affy probe ids, and a csv file from affymetrix with
    the appropriate information
    @param affy_id_list: A list of strings like '1000_at'...
    @param affy_file: If none, then the function get_affy_key_file()
    will be called to get the full file name and path to the csv file,
    else specify the filename/path.
    '''
    if affy_file is None:
        affy_file = get_affy_key_file()
        
    affy_dict = {}
    
    lines = []
    with open(affy_file,"r") as f:
        for tmpline in f:
            if tmpline[0] != "#": lines.append(tmpline)  #omit header/comment lines
            
    for i,ln in enumerate( lines[1:] ):  #lines[0] is the column headers
        ifr.print_progress(i, len(lines))
        tmp = ifr.smart_split(ln, sep=",")
        key = tmp[0]
        genebank = tmp[8]
        unigene = tmp[10]
        symbol = tmp[14]
        affy_dict[ key ] = (genebank,unigene,symbol)
        
    return affy_dict

def fix_affy_dict(affy_dict):
    '''
    Given an affy dict in form { probe_id:(genebank,unigene,symbol)},
    such as generated by get_geneids_from_affy() by reading a CSV file,
    this attempts to fix those entries where the symbol was missing or
    multiple-valued by doing a live query against the latest symbol
    information from ncbi web sources.
    '''
    #N = len(affy_dict)
    new_dict = {}
    for ix,k in enumerate( affy_dict.keys() ):
        lookup = False
        #svo_util.print_progress(ix, N)
        (genebank, unigene,symbol) = affy_dict[k]
        if symbol == '---':
            print "(%d) %s unknown"%((ix+1),k) , 
            lookup = True
        elif '///' in symbol:
            print "(%d) %s ambiguous"%((ix+1),k) , 
            lookup = True
            
        if lookup:
            if (unigene == '---') or ('///' in unigene):
                symbol = ifr.get_official_name(genebank)
            else:
                #if we have a unigene name, try it first,
                # but if it doesn't work, then try the genebank
                symbol = ifr.get_official_name(unigene)
                if symbol == unigene:
                    symbol = ifr.get_official_name(genebank)
            print " = %s"%symbol
            
        new_dict[k] = (genebank, unigene, symbol)
    
    return new_dict

def get_genelist_from_affy( affy_id_list, d=None):
    '''
    Returns a list of gene ids given an input list of affy probe ids.
    @param affy_id_list: The input list of affy probe ids, like ['1000_at', ...]
    @param d: The mapping dictionary, if already computed. This could be the
    result of the functions load_affy_to_geneId_dict() or gen_affy_to_geneId_dict().
    If None, then the dict will be loaded from the default pickle file location.
    @return: A list of gene ids like ['MAPK3', ...]    
    @note: There is not a 1-to-1 relationship. Some affy probe ids may have multiple
    gene symbols listed. Thus the length of the returned list may not be the same
    as the input list.
    '''
    affy_id_list = [ aid.strip() for aid in affy_id_list]  #remove leading/trailing spaces
    
    if d is None:
        d = load_affy_to_geneId_dict()
        
    genelist = []
    for a in affy_id_list:
        if a in d:
            genelist += d[a]
        else:
            print "Warning: affy_id %s not found in dictionary."%a
        
    return sorted( list(set(genelist)) ) #remove duplicates
    

def get_affy_key_file():
    '''
    Return the file name of the affymetrix csv file that defines the probe sets.
    This is just a convenience function, assumes the file subdir and name are default.
    '''
    return os.path.join( ifr.DATA_DIR, "HG_U95A.na33.annot", "HG_U95A.na33.annot.csv")
    
def load_affy_to_geneId_dict(affy_file_subdir="HG_U95A.na33.annot", dict_fn="affy_to_gene_id_dict.p"):
    '''
    Convenience function for loading the affy probe id to gene id dictionary from a saved pickle file.
    The dictionary should be as generated by gen_affy_to_geneId_dict() function, but it is more
    efficient to save the dictionary than re-parse the csv file.
    '''
    fn = os.path.join(ifr.DATA_DIR, affy_file_subdir, dict_fn)
    return cPickle.load( open(fn,"rb") )
    
def gen_affy_to_geneId_dict(affy_file_subdir="HG_U95A.na33.annot", affy_fn="HG_U95A.na33.annot.csv"):
    '''
    Converts a list of affymetric probe set ids into a genelist
    with names suitable for querying gather or kegg. Generates
    a dictionary with entries { affy_id : gene_id_list }. Most times, gene_id_list will
    have only a single entry, but several probes have multiple Gene IDs given.
    @param affy_file_subdir: The subdirectory of the ifr.DATA_DIR that has
    the HG_U95A.na33.annot.csv file.
    @param affy_fn: The csv file in the subdirectory with the data. The parameter is
    provided in case the file was renamed from the orginal name of "HG_U95A.na33.annot.csv"
    @note: Relies on a data file called HG_U95A.na33.annot.csv that must be present
    in the HG_U95A.na33.annot subdirectory of the linked Data directory. It would
    be most efficient to use this function once and save the resulting dictionary
    in a pickle file for later use instead of having to re-parse the data.
    '''
    affy_dict = {}
    
    affy_file = os.path.join( ifr.DATA_DIR, affy_file_subdir, affy_fn)
    lines = []
    with open(affy_file,"r") as f:
        for tmpline in f:
            if tmpline[0] != "#": lines.append(tmpline)  #omit header/comment lines
            
    for i,ln in enumerate( lines[1:] ):  #lines[0] is the column headers
        ifr.print_progress(i, len(lines))
        tmp = ifr.smart_split(ln, sep=",")
        key = tmp[0]
        val = tmp[14]
        if val == '---':
            #this affy id has no gene symbol
            genelist = []
        elif "///" in val: #there are more than one GeneIds for this probe
            #print "Subfield indicator in Gene Symbol for %s, line: %d."%(key,(i+1))
            #print "Val: %s"%tmp[14]
            genelist = [ x.strip() for x in val.split("///") if x.strip() != '' ]
        else:
            genelist = [val]
            
        affy_dict[ key ] = list(set(genelist)) #remove duplicates
        
    return affy_dict
