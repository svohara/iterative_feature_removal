'''
Created on Oct 29, 2012
@author: ohara

Script interface for querying Duke's GATHER system
to present annotations from the gene ontology for a
set of genes

'''
import ifr
import math

class gather(object):
    def __init__(self, url='http://gather.genome.duke.edu'):
        self.url = url
        self.data = []
        
    def query(self, genelist, homologs=True):
        '''
        Processes the query for the given list of genes. The result will be to cache the
        results in self.data. Use other methods of this class to query the cached data, for
        example, to get the data for annotation GO:0009615, which is 'response to virus',
        use getDataForAnnotation()
        @param genelist: A list of gene strings like ['RSAD2', 'SERPING1',....]
        @param homologs: If True, then the query includes annotations for homologs to the input
        genelist, else it does not.
        @return: A return code indicating status of request. The data will be
        stored in self.data[] as a list of lines.
        '''
        self.data = []
        gene_box = ",".join(genelist)
        h = '1' if homologs else '0'
        p = {'annot_type': 'gene_ontology', 'cmd': 'report', 'homologs': h, 'network': '0', 'tax_id': '9606'}
        p['gene_box'] = gene_box
        
        raw_data, rc = ifr.fetch_data(self.url, p, method='POST')  #post doesn't work?
        for ln in raw_data:
            self.data.append( ln.split('\n')[0])
                
        return rc
    
    def getDataForAnnotation(self, key='GO:0009615', return_mult=True, bf_thresh=0.0):
        '''
        Returns the parsed data for all rows which contain a given key,
        which can be any string that one expects to find in the response data from a gather query.
        @param key: A string, which could be a gene ontology id, a word (such as 'virus'), or any
        other text that can be found in the result set.
        @param return_mult: Set true to return multiple rows which meet the search key. False
        finds only the first row that has the search key.
        @param bf_thresh: The bayes factor threshold so that only annotations meeting or exceeding
        this threshold will be returned.
        @return: A tuple (fields, dat), where fields provides the header data defining the fields
        of the returned data, and dat is a list of lists. Each top-level list is from a row of the response
        that contains the key string. The second level is the list of field values for each row, where
        the field values correspond to the 'fields' header information.
        @note: You must have called the query() method to load the gather data for
        the query genes prior to using this function.
        '''
        if len(self.data) < 1:
            print "Error: Data has not been loaded. Call the query() method first."
            return
        
        dat = []
        for ln in self.data[1:]: #first line is header titles
            line_dat = ln.split('\t')
            bf_tmp = float(line_dat[7])
            if key in ln:
                #print bf_tmp, bf_thresh, bf_tmp >= bf_thresh
                if bf_tmp >= bf_thresh:
                    dat.append(ln.split('\t'))
                    if not return_mult: break  #do we expect only a single result/key? if so, we're done searching
                    
        if len(dat) < 1:
            print "Error: Requested key: %s is not in result data with high enough bayes factor score."%key
            return
        
        #first line in self.data should be the field definitions
        fields = self.data[0].split('\t')
                
        return (fields,dat)
    
    def getRowIdsForAnnotation(self, key='virus'):
        '''
        Returns all the row ids for the current query results where the data
        contains the string provided by key.
        '''
        if not type(key)==list:
            key = [key]
            
        if len(self.data) < 1:
            print "Error: Data has not been loaded. Call the query() method first."
            return
        
        dat = []
        for ln in self.data:
            for k in key:
                if k in ln:
                    rowId = int(ln.split('\t')[0])
                    if not rowId in dat: dat.append(rowId)
            
        if len(dat) < 1:
            print "Error: Requested key: %s is not in result data."%key
            return
                
        return dat
        
    
    def getBayesFactor(self, row_id):
        '''
        Given the gather row identifier (not the index of self.data, but rather the first field
        of a response row in the data), returns the bayes factor and p-value.
        @param row_id: The row identifier, which should be an integer.
        @return: (label, bayes factor, p-value). See the Gather documentation for an explanation.
        @note: self.data must be populated by a call to query() before this method will work.
        '''
        if len(self.data) < 1:
            print "Error: Data has not been loaded. Call the query() method first."
            return
        
        bf_tmp = None
        
        for line in self.data[1:]: #first line is header info
            line_dat = line.split('\t')
            line_id = line_dat[0]
            if int(line_id) == row_id:
                bf_tmp = line_dat[7]
                pv_tmp = line_dat[8]
                label = line_dat[1]
                break
        
        if bf_tmp is None:
            print "Error: There is no row #: %d."%row_id
            return (None, None)
        
        bf = float(bf_tmp)
        pv = 1.0 / math.exp( float(pv_tmp))
        
        return (label, bf, pv)
    
    def parseAnnotation(self, annotationList):
        '''
        Parses an annotation to generate a dictionary
        with key-value bindings.
        @return: a dictionary with key/value bindings describing the fields
        of a given annotation.
        '''
        d = {}
        d['01. RowId'] = int(annotationList[0])
        tmp = annotationList[1].split(' ')
        
        d['02. AnnotationId'] = tmp[0]
        d['03. OntologyDepth'] = int(tmp[1][1:-2])
        d['04. Description'] = (" ".join( tmp[2:] )).strip()
        
        d['05. Total Genes With Ann'] = int( annotationList[2])
        d['06. Your Genes With Ann'] = int( annotationList[3])
        d['07. Your Genes No Ann'] = int( annotationList[4])
        d['08. Genome With Ann'] = int( annotationList[5])
        d['09. Genome No Ann'] = int( annotationList[6])
        d['10. ln(Bayes Factor)'] = float( annotationList[7])
        d['11. neg ln(p val)'] = float( annotationList[8])
        d['12. FE: neg ln(p val)'] = float( annotationList[9])
        d['13. FE: neg ln(FDR)'] = float(annotationList[10])
        d['14. Genes'] = annotationList[11]
        
        return d
    
    def countRows(self, bf_thresh=0):
        '''
        counts the number of data rows (annotations)
        that have a bayes factor score >= bf_thresh
        '''
        if len(self.data) < 1:
            print "Error: Data has not been loaded. Call the query() method first."
            return
        
        count = 0
        for line in self.data[1:]:
            line_dat = line.split('\t')
            if len(line_dat) < 8: continue  #bad line?
            bf = float(line_dat[7])
            if bf >= bf_thresh: count += 1
            
        return count

def test_gather(key="virus"):
    genelist = ['LY6E','SERPING1','IFI44L','IFI27','OAS1']
    g = gather()
    g.query(genelist)
    (_fields, dat) = g.getDataForAnnotation(key=key)
    return dat
        
if __name__ == '__main__':
    pass