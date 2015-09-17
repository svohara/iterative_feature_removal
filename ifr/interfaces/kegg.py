'''
Created on Dec 21, 2012
@author: Stephen O'Hara

Interface for querying KEGG web service data
given gene names
'''
from ifr import fetch_data

def kegg_get_gene_id(g="CDK5", species_code="hsa"):
    '''
    Return the KEGG gene identifier for a given common gene symbol.
    @param g: The gene common name, like "OAS1"
    @param species_code: The code for the genome where this gene is from. See
    KEGG for more details. hsa is the human genome, mmu is mouse.
    @return: The kegg gene id.
    '''
    rc = fetch_data("http://rest.kegg.jp/find/%s/%s"%(species_code,g), {}, method='POST')
    if rc[1] != 200:
        print "Query failed with error code: %d"%rc[1]
        return None
    
    if len(rc[0]) < 1 or rc[0] == ['\n']:
        print "Query succeeded, but has no data for gene: %s"%g
        return None
    
    kegg_id = None
    for line in rc[0]:                                     
        tmp = line.split('\t')[1].split(';')[0].split(',')
        tmp = [ t.strip() for t in tmp]
        #print "Debug: ", tmp
        if g in tmp:
            kegg_id = line.split('\t')[0].strip()
            break
        
    return kegg_id
        
def kegg_get_pathway_info(p="path:hsa05222"):
    '''
    Returns all the information for a given kegg pathway identifier
    '''
    rc = fetch_data("http://rest.kegg.jp/get/%s"%p, {}, method='POST')
    if rc[1] != 200:
        print "Query failed with error code: %d"%rc[1]
        return None
    
    if len(rc[0]) < 1 or rc[0] == ['\n']:
        print "Query succeeded, but has no data for pathway: %s"%p
        return None
    
    return rc[0]
    
def kegg_get_genes_from_pathway(p="path:hsa05222"):
    '''
    Returns a list of gene symbols that are listed as being
    included in a given kegg pathway
    '''
    info = kegg_get_pathway_info(p)
    if info is None:
        return None
    
    #find where the gene list starts
    start_idx = 9999
    for i,line in enumerate(info):
        if line[0:4] == 'GENE':
            start_idx = i 
            break
    assert start_idx < 9999 
    
    #find where gene list ends
    end_idx = -1
    for i,line in enumerate( info[(start_idx+1):] ):
        if line[0:12] != (' '*12): #twelve spaces leader
            end_idx = i + start_idx
            break        
    assert end_idx >= start_idx
    
    genelines = info[start_idx:(end_idx+1)]
    genelines[0] = genelines[0][4:]  #get rid of the GENE tag on first line
    genelines = [ gl.strip() for gl in genelines]  #get rid of useless whitespace in all rows
    
    genelines2 = [ x.split(' ')[2].split(';')[0] for x in genelines] #just the gene name parsed from the rest
    
    return genelines2
            
def kegg_get_pathway_name(p="path:hsa05164"):
    '''
    Returns the kegg pathway name (description) for a given path identifier
    @param p: The kegg pathway identifier
    '''
    rc = fetch_data("http://rest.kegg.jp/list/%s"%p, {}, method='POST')
    if rc[1] != 200:
        print "Query failed with error code: %d"%rc[1]
        return None
    
    if len(rc[0]) < 1 or rc[0] == ['\n']:
        print "Query succeeded, but has no data for pathway: %s"%p
        return None
    
    return rc[0][0].split('\t')[1].strip()
        
def kegg_get_pathways_with_gene(g="OAS1", species_code="hsa"):
    '''
    Determines the set of known pathways associated to a specific gene
    @param g: The common gene symbol
    @param species_code: Like "hsa" for human, "mmu" for mouse. See the
    KEGG documentation for available organism codes.
    @return: Tuple (kegg_path_names, kegg_path_ids)
    '''
    kegg_id = kegg_get_gene_id(g, species_code)
    if kegg_id is None:
        print "Error: Could not resolve KEGG path id for gene name %s"%g
        return None
    
    rc = fetch_data("http://rest.kegg.jp/link/pathway/%s"%kegg_id, {}, method='POST')
    if rc[1] != 200:
        print "Query failed with error code: %d"%rc[1]
        return None
    
    if len(rc[0]) < 1 or rc[0] == ['\n']:
        print "Query succeeded, but has no pathway ids that include gene: %s"%g
        return None
    
    kegg_path_ids = [ line.split('\t')[1].strip() for line in rc[0] ]
    kegg_path_names = []
    for pid in kegg_path_ids:
        tmp = kegg_get_pathway_name(pid)
        if tmp is None: tmp = "unknown"
        kegg_path_names.append(tmp)
    
    return (kegg_path_names, kegg_path_ids)
        
def kegg_get_pathways_with_genelist(genelist, species_code="hsa"):
    '''
    Determines a set of known pathways associated with a list of genes. This
    is done by repeatedly calling kegg_get_pathways_with_gene() function for
    each gene in the list, and combining results where a pathway is listed for
    more than one of the genes in the list.
    @param genelist: A list of gene names
    @param species_code: The KEGG organism code to search, "hsa" is human, "mmu" is mouse.
    @return: A dictionary object, where the keys are pathway ids, and the values are tuples
    of format (description, num_genes). Description is the descriptive label for the pathway,
    and num_genes indicates how many of the genes in genelist were associated to this pathway.
    '''
    pathway_dict = {}  # {pathway_id:[desc, num_genes]}
    
    for g in genelist:
        rc = kegg_get_pathways_with_gene(g, species_code)
        if rc is None: continue
        
        for (d,pid) in zip(*rc):
            if pid in pathway_dict:  #already known for this genelist, but increment count
                tmp = pathway_dict[pid]
                tmp[1] += 1
            else:              
                pathway_dict[pid] = [d,1]
    
    return pathway_dict
    
def compute_kegg_annotations(genelist,species_code="hsa"):
    '''
    Computes the data required to produce a bar chart for the genelist where the y-axis has
    kegg pathway annotations, and the x-axis is the count of genes in genelist that have that annotation.
    '''
    pathway_dict = kegg_get_pathways_with_genelist(genelist,species_code)
    if len(pathway_dict) < 1:
        print "Error: No pathway annotations for given gene list."
        return []
    
    #dat = []
    #for pid in pathway_dict:
    #    [desc, count] = pathway_dict[pid]
    #    dat.append( (pid, desc, count))
    
    dat = sorted( [ (pid,desc,count) for (pid,(desc,count)) in pathway_dict.items()] )
    return dat 
    
        
if __name__ == '__main__':
    pass