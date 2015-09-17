'''
Created on Dec 27, 2012
@author: Stephen O'Hara

Convenience functions for loading/reading certain source data files.
'''
import os
import scipy.io as spio
import scipy as sp
import cPickle
import ifr
import zipfile

class CSV_data:
    '''
    Class to encapsulate reading data from simple csv files
    where the first row is the column headers and fields are
    separated with commas. Rows are data samples, columns are field values.
    '''
    def __init__(self, datadir, csvfile, label_cols=[0], separator="," ):
        '''
        constructor
        @param datadir: The directory where the csv file is located
        @param csvfile: The filename of the data file
        @param label_cols: A list of the zero-indexed column numbers
        that should be treated as labels. These columns are assumed
        to have a discrete set of values, like a set of string labels,
        a set of integers. For regression data sets where the "labels"
        are continuous-valued numbers, set label_cols = None, and all
        the csv fields will be loaded into the data matrix.
        '''
        self.datadir = datadir
        self.csvfile = csvfile
        self.label_cols = label_cols
        self.label_dict = {}  # key is label_col #, value is a list of the values of each sample
        self.label_names = [] #list of the field names of the label columns
        self.fields = [] #list of the field names for the data columns (non-labels)
        self.separator = separator
        self.skipped_rows = []  #row skipped if a data field was empty
        
        self._load_data()
        
        
    def _load_data(self):
        '''
        internal function that reads the csv text data and converts
        the values as appropriate to the internal data representation
        '''
        infile = os.path.join(self.datadir,self.csvfile)
        with open(infile,"r") as f:
            lines = f.readlines()
        first_row = self._load_header(lines)
        raw_dat = lines[(first_row+1):]  #drop the header row and any leading blank rows
        
        data_list = []   
        for i,linedat in enumerate(raw_dat):
            row_num = i + first_row + 1 #original row from csv file including header/initial blank lines
            row_dat = self._parse_row(linedat)
            if row_dat is None:
                self.skipped_rows.append(row_num)
                print "Warning: row %d of input skipped with missing values."%row_num
                #print "Line data: %s"%str(linedat)
            else:
                data_list.append(row_dat)
         
        self.data = sp.array(data_list)
        print "Loaded data matrix of size: %d by %d"%self.data.shape
        print "The header row was line %d"%(first_row+1) #one-based for user
        print "The label columns are: %s"%str(self.label_names)
        
    def _load_header(self, lines):
        '''
        header has the field names, and should
        be the first non-empty line of the csv file
        '''
        i = 0
        while lines[i].strip() == '':
            i+=1
        first_row = i
        self.fields = [ fld.strip() for fld in lines[first_row].split(self.separator) ] 
        
        #setup initial label_dict entries for label fields
        for col in self.label_cols:
            self.label_dict[ col ] = []
            self.label_names.append(self.fields[col])
            
        #self.fields will have the field names for the data (non-label) columns in order
        for label_name in self.label_names:
            self.fields.remove(label_name)
                
        return first_row
    
    def _parse_row(self, linedat):
        linedat = linedat.strip()
        tmpfields = [ fld.strip() for fld in linedat.split(self.separator)]
        
        #are any field values empty? then badrow        
        if '' in tmpfields: return None
        
        #handle label columns
        for col in self.label_cols:
            fieldval = tmpfields[col]
            self.label_dict[col].append( fieldval )
            
        #non-label columns are converted to floats, appended into an ordered list
        row_dat = []
        data_cols = set(range(len(tmpfields)))-set(self.label_cols)
        for x in sorted(list(data_cols)):
            row_dat.append( float(tmpfields[x]) )
            
        return row_dat

class C45_data:
    '''
    Class to encapsulate functions to read data from C4.5 formatted
    data mining files. Characterized by text files with .names and .data
    extensions.
    '''
    def __init__(self, datadir, namefile, datafile, missing_val_sentinel=-9999):
        '''
        Constructor
        @param datadir: The directory where the data files are located
        @param namefile: The filename that has the .names information. You may
        specify as xyz.names or just xyz (with .names extension inferred)
        @param datafile: The filename that has the .data information. You may
        specify as xyz.data or just xyz (with .data extension inferred)
        @param missing_val_sentinel: Set to a numeric value that can uniquely identify
        fields in the data matrix where the values were missing. Matrix can be
        transformed via functions like 'replace_missing_values_with_col_means()'. Set
        to None to have an error thrown if any missing values are encountered.
        @note: Do NOT include the path in the namefile or datafile parameters. These
        are assumed to be files in the datadir specified.
        '''
        self.classes = []
        self.fields = []
        self.fieldtypes = []
        self.data = []
        self.labels = []
        
        self.datadir = datadir
        self.namefile = namefile if namefile[-5:] == "names" else "%s.names"%namefile
        self.datafile = datafile if datafile[-4:] == "data" else "%s.data"%datafile
        
        self.missing_val_sentinel = missing_val_sentinel
        
        self._loadNames()
        self._loadData()
        
    def _loadNames(self):
        nfile = os.path.join(self.datadir, self.namefile)
        assert os.path.exists(nfile)
        
        self.classes = []
        self.fields = []
        self.fieldtypes = []
        
        with open(nfile,"r") as f:
            lines = f.readlines()
            
        #first line is supposed to be the class names
        tmp = lines[0].split(",")
        self.classes = [c.strip() for c in tmp]
        print "There are %d classes defined by data set: %s"%(len(self.classes), str(self.classes))
        
        for tmp2 in lines[1:]:
            #if the line is just whitespace, skip it
            if tmp2.strip() == '': continue
            tmp3 = tmp2.split(':')
            self.fields.append( tmp3[0].strip())
            self.fieldtypes.append( tmp3[1].strip() )
        
        print "There are %d fields in the data set."%len(self.fields)
        
    def _parse_row(self, tmpfields):
        '''
        Given a list of strings representing one row of data, return a vector of numbers.
        This will replace missing values with the sentinel number, or raise an error if
        self.missing_val_sentinel is None.
        '''        
        #are any field values 'empty' or non-numeric
        try:
            V = sp.array(tmpfields, dtype=float)
        except(ValueError):
            #There are values in tmpfields that can't be converted
            # to floats, we assume these are missing values.
            if self.missing_val_sentinel is None:
                raise(ValueError)
            else:
                #replace any non-numeric entries in tmpfields with the sentinel number
                tmp2 = [ ifr.parse_number(s, fail=self.missing_val_sentinel) for s in tmpfields ]
                V = sp.array(tmp2, dtype=float)
                    
        return V
       
    def _loadData(self):
        '''
        Loads the data. self.data will be a numpy matrix representing the source data.
        Each row is a sample, and each column is a feature.
        '''
        dfile = os.path.join(self.datadir, self.datafile)
        assert os.path.exists(dfile)
        data = []
        labels = []
        
        with open(dfile, "r") as f:
            datalines = f.readlines()
            
        for dl in datalines:
            tmp = dl.split(',')
            V = self._parse_row( tmp[:-1] )
            data.append( V )
            labels.append(tmp[-1].strip())  #last entry for each row is the label
            
        print "Loaded %d lines of data."%len(data)
        self.data = sp.array(data, dtype=float)
        self.labels = sp.array(labels)
        P = self.data.shape[1]
            
        #some of the data features may be constant (stddev=0), so they need to be dropped
        cc = ifr.constant_cols(self.data)
        if len(cc) > 1:
            print "Warning, %d features in data set have zero standard deviation."%len(cc)
            print "Those features will be dropped."
            all_cols = set(range(P))
            good_cols = sorted( list(all_cols - set(cc))) #set difference
            #filter field names to remove bad columns
            self.fields = [ f for i,f in enumerate(self.fields) if i in good_cols]
            self.fieldtypes = [ f for i,f in enumerate(self.fieldtypes) if i in good_cols]
            #filter data array to remove bad columns
            mask = sp.array(good_cols)
            self.data = self.data[:,mask]
            #remember which were the bad columns as an object property
            self.constant_cols = cc
    
    def labelsAsIntegers(self, include_label_set=False):
        '''
        Returns the labels as a list of integer values, and the accompanying ordered
        list that tells you which text label is represented by which integer
        @param include_label_set: If True, then the return will be a tuple of
        the label integers and a list of the unique labels in the index order
        of the integers. If False, just the label integers are provided in a single
        value return.
        @return: either int_list or tuple (int_list, unique_labels)
        '''
        unique_labels = sorted(list(set(self.labels)))
        integer_labels = [ unique_labels.index(lbl) for lbl in self.labels]
        if include_label_set:
            return (integer_labels , unique_labels)
        else:
            return integer_labels
        
    def labelsAsIndicatorVals(self, include_label_set=False):
        '''
        Returns the labels as indicator variables. Indicator variables have a column for
        each possible class, and either a 0 or 1 in that column if that label applies.
        It is common for classification algorithms to take indicator variables as input
        for training data labels.
        @return: tuple ( IV_matrix, unique_labels )
        '''
        RC = ifr.indicator_vars(self.labels)
        if include_label_set:
            return RC
        else:
            return RC[0]


def install_kentridge_data(data_set):
    '''
    Downloads and installs (unzips/extracts) one of the kent ridge data sets.
    @param data_set: Use ifr.LUNG_CANCER_URL, ifr.PROSTATE_URL, or ifr.BCELL_LYMPHOMA_URL
    @return: This function has no return. When finished without error, the specified
    kent ridge data should be downloaded and extracted and ready to use.
    '''
    allowed_urls = (ifr.LUNG_CANCER_URL, ifr.PROSTATE_URL, ifr.BCELL_LYMPHOMA_URL)
    assert data_set in allowed_urls, "Error: you must specify one of the three known data sets."

    #Note in below...the prostate zip file is different than the other two because it has
    # the data in a subdir called prostate, so the top dir should be the kent_ridge dir.
    output_paths = (ifr.LUNG_CANCER_DIR, ifr.KENT_RIDGE_DIR, ifr.BCELL_LYMPHOMA_DIR)
    out_path = output_paths[ allowed_urls.index(data_set)]

    zipfilename = os.path.basename(data_set)
    dest_file = os.path.join(ifr.KENT_RIDGE_DIR, zipfilename)
    
    print "Downloading zip archive from %s"%data_set
    rc = ifr.download_file(data_set, dest_file)
    
    print 'Extracting data from zip archive'
    z = zipfile.ZipFile(rc[0],'r')
    z.extractall(out_path)
    
    return rc
    

def load_flu_mat(matfile="Fluz_dat.mat", times=[14,15,16], id_filter=None):
    '''
    Loads the Duke flu data, only those rows where the response
    is for the given time indexes.
    @param matfile: A matlab file that contains duke influenza data,
    but converted to a format easily read by python/scipy. The .mat file
    should contain four matrices (not cell arrays), D, L, Time, ID,
    where D is the full data, L is the corresponding labels, Time is the time
    index of the data, and ID is the subject id for each row of data.
    @param times: A list of time indexes that is used to filter the
    data. Only those rows corresponding to these time indexes will be
    included. Specify None to include all time steps.
    @param id_filter: A list of ids to FILTER (remove) from the data.
    Specify None to include all subjects.
    @return: tuple (D,L,Time,ID) as numpy arrays with the filtered data.
    @note: For Duke Fluz data, use times=[14,15,16] and id_filter=None. For
    Duke H1N1 data, use times=[14,15,16] and id_filter = [1,3,5,10,18]
    '''
    print "Loading data from %s in directory: %s"%(matfile, ifr.DUKE_PYDAT_DIR)
    fn = os.path.join(ifr.DUKE_PYDAT_DIR, matfile)
    tmp = spio.loadmat(fn)
    D = tmp['D']  #full data
    L = tmp['L'].flatten()  #all labels
    Time = tmp['Time'].flatten()  #indicates the time index for each row
    Ids = tmp['ID'].flatten() #subject id for each row
    if not (times is None):
        idxs = []
        for t in times:
            idxs += list(sp.flatnonzero(Time==t))
        
        D = D[idxs,:]
        L = L[idxs,:]
        Time = Time[idxs,:]
        Ids = Ids[idxs,:]
        
    if not (id_filter is None):
        idxs = []
        for i in list(sp.unique(Ids)):
            if not (i in id_filter):
                idxs += list(sp.flatnonzero(Ids==i))
        
        D = D[idxs,:]
        L = L[idxs,:]
        Time = Time[idxs,:]
        Ids = Ids[idxs,:]
    
    return (D,L,Time,Ids)        

def load_H1N1_mat():
    '''
    convenience function to call load_flu_mat() with arguments set
    to load the 57 samples of the Duke H1N1 data
    '''
    return load_flu_mat('H1N1_dat.mat', times=[14,15,16], id_filter=[1,3,5,10,18])

def load_housekeeping_genes(G):
    '''
    G is a list of genes from the Duke data set (affymetrix).
    This function returns a list of nominal housekeeping genes,
    which are those starting with a particular prefix.
    '''
    #perhaps all genes starting 'AFFX-....' are housekeeping,
    # but this is a set of 21 that seem to work well
    hk_genes = [ g for g in G if g[0:7]=='AFFX-r2']
    return hk_genes

def load_gene_ids(fn="gene_ids.csv", ddir=ifr.DUKE_PYDAT_DIR, short_name_only=False):    
    '''
    Loads the gene names/identifiers from a csv file
    '''
    csvf = os.path.join(ddir,fn)
    with open(csvf,"r") as f:
        dat = f.readlines()
        
    gene_ids = []
    for line in dat:
        x = line.split(",")
        id1 = x[0]
        id2 = x[1].split("\n")[0]
        gene_ids.append((id1,id2)) 
    
    if short_name_only:
        gene_ids = [ gid for (_,gid) in gene_ids]
    
    return gene_ids

def load_lung_cancer_dat():
    '''
    Loads the lung cancer data downloaded from the kent ridge biomedical repository.
    The data is in C45 format, and is comprised of two files, the train and test data.
    @return: (train_c45_object, test_c45_object).
    '''
    if not os.path.exists( os.path.join(ifr.LUNG_CANCER_DIR, "lungCancer_train.data")):
        print "This appears to be the first call to loading the lung cancer data."
        print "The data will now be downloaded and installed. An internet connection is assumed..."
        _rc = install_kentridge_data(ifr.LUNG_CANCER_URL)
        
    x_train = ifr.C45_data(ifr.LUNG_CANCER_DIR, "lungCancer_train", "lungCancer_train")
    x_test = ifr.C45_data(ifr.LUNG_CANCER_DIR, "lungCancer_test", "lungCancer_test")
    return (x_train, x_test)

def load_lung_cancer_dict():
    '''
    Loads a pre-computed dictionary that maps the lung cancer's
    affy probe ids to other gene ids (genebank,unigene,symbol)
    that are more useful.
    '''
    fn = os.path.join('ifr','interfaces','affyU95_dict_lungcancer.p')
    return cPickle.load(open(fn,'rb'))

def load_prostate_dat():
    '''
    as per load_lung_cancer_dat, same format but different data set
    @note: Loads the "Tumor vs Normal" section of the data set.
    '''
    tmp = "prostate_TumorVSNormal"
    
    if not os.path.exists( os.path.join(ifr.PROSTATE_DIR, "%s_train.data"%tmp)):
        print "This appears to be the first call to loading the prostate cancer data."
        print "The data will now be downloaded and installed. An internet connection is assumed..."
        _rc = install_kentridge_data(ifr.PROSTATE_URL)
    
    x_train = ifr.C45_data(ifr.PROSTATE_DIR, tmp, "%s_train"%tmp)
    x_test = ifr.C45_data(ifr.PROSTATE_DIR, tmp, "%s_test"%tmp)
    return (x_train, x_test)

def load_bcell_lymphoma_dat():
    '''
    as per load_lung_cancer_dat, same format but different data set
    @note: has missing values, which are replaced by this code with
    the column means. The field names are like GENE2148X, where 2148
    is the entrezgene identifier, which can be converted to a gene
    symbol using ifr.get_official_name( '2148', entrezgene=True)
    '''
    if not os.path.exists( os.path.join(ifr.BCELL_LYMPHOMA_DIR, "DLBCL.data")):
        print "This appears to be the first call to loading the bcell lymphoma data."
        print "The data will now be downloaded and installed. An internet connection is assumed..."
        _rc = install_kentridge_data(ifr.BCELL_LYMPHOMA_URL)
        
    x = ifr.C45_data(ifr.BCELL_LYMPHOMA_DIR, "DLBCL", "DLBCL")
    D = ifr.replace_missing_values_with_col_means(x.data, sentinel=-9999)
    x.data = D
    return x

def load_bcell_lymphoma_dict():
    '''
    Loads a dictionary that maps the field names from the
    bcell lymphoma data set to gene symbols
    '''
    fn = os.path.join('ifr','interfaces','bcell_fields_dict.p')
    return cPickle.load(open(fn,'rb'))


    
    
    
    
    
    
