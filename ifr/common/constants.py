'''
Created on Sep 25, 2012
@author: Stephen O'Hara

Handy constants
'''
import os
import data as ifrdata

DATA_DIR = ifrdata.__path__[0]
RESULTS_DIR = os.path.join(DATA_DIR, "results")

#directory where the python-compatible data files are located
# created from the Duke data
DUKE_PYDAT_DIR = os.path.join(DATA_DIR, "Duke_Data_For_Python")

#DUKE Panviral datasets
PV_DATASET_H1N1 = 'H1N1'
PV_DATASET_H3N2 = 'Fluz'
PV_DATASET_HRV = 'Rhino'
PV_DATASET_RSV = 'RSV'
PV_DATASETS = [PV_DATASET_H3N2, PV_DATASET_H1N1, PV_DATASET_HRV, PV_DATASET_RSV]

#Data sets from Kent Ridge repository
KENT_RIDGE_DIR = os.path.join(DATA_DIR, "KentRidgeData")
LUNG_CANCER_DIR = os.path.join( KENT_RIDGE_DIR, "lungcancer")
PROSTATE_DIR = os.path.join( KENT_RIDGE_DIR, "prostate")
BCELL_LYMPHOMA_DIR = os.path.join( KENT_RIDGE_DIR, "DLBCL")
LUNG_CANCER_URL = 'http://levis.tongji.edu.cn/gzli/data/lungcancer.zip'
PROSTATE_URL = 'http://levis.tongji.edu.cn/gzli/data/Prostate.zip'
BCELL_LYMPHOMA_URL = 'http://levis.tongji.edu.cn/gzli/data/DLBCL.zip'


#Files output by matlab SSVM-based IFR
IFR_INFLUENZA_FILE_GENES = os.path.join(RESULTS_DIR, "IFR_SSVM", "H3N2_Genelist_SSVM_02112013.csv")
IFR_INFLUENZA_FILE_ACC = os.path.join(RESULTS_DIR, "IFR_SSVM", "H3N2_SSVM_Remove_Acc_02112013.txt")
IFR_LUNG_FILE_GENES = os.path.join(RESULTS_DIR, "IFR_SSVM", "Lung_Genelist_SSVM_02112013.csv")
IFR_LUNG_FILE_ACC = os.path.join(RESULTS_DIR, "IFR_SSVM", "Lung_SSVM_Remove_Acc_02112013.txt")

class flu_genelists:
    '''
    Structure for organizing all the constant genelists that are used to generate various
    figures.
    '''
   
    #The following list from Chen et al, 2011, page 8. Note that 'AFFX-HUMISGF3A/M97935_3_at' was reported as HU-MISGF3A
    BENET50_GENELIST = ['OAS1', 'IFI44', 'DDX60', 'IFI44L',  'IFIT1', 'GBP1', 'LY6E', 'RSAD2', 'SCO2',
                             'TRIM22', 'IFITM1', 'IFIT3', 'ISG15', 'MS4A4A', 'IFI6', 'PSME2', 'XAF1', 'IFIT5', 'IFITM3',
                             'HERC5', 'IFI27', 'SERPING1', 'PLSCR1', 'OAS3', 'STAT1', 'ZCCHC2', 'AFFX-HUMISGF3A/M97935_3_at', 'OASL',
                             'C1QB', 'SIGLEC1', 'AIM2', 'OAS2', 'MX1', 'LOC26010', 'SMPDL3A', 'C13orf18', 'LAP3', 'TDRD7',
                             'PARP12', 'PSME1', 'VAMP5', 'IRF7', 'SAT1', 'SMAD1', 'BLVRA', 'IDH2', 'C1QA', 'MT2A', 'TAP1']
    
    
    #Interferon Stimulated Genes found in the first 40 iterations of SSVM IFR
    # The AFFX-... are described as relating to the STAT1 pathway and thus should? be part of this group.
    SSVM_ISG = ['OAS1','OAS3','OASL', 'OAS2', 'RSAD2', 'MX1', 'STAT1', 'STAT2', 'STAT5B', 'HERC5', 'HERC6',
                'IFI44', 'IFIT1', 'IFI44L', 'IFITM1', 'IFI27', 'IFIT3', 'IFI6', 'IFITM3', 'IFIT5', 'IFIH1', 'IFI35',
                'IFI30', 'IFIT2', 'IFITM2', 'ISG15', 'IL15', 'IL15RA',
                'DDX60', 'GBP1', 'TRIM22', 'PSME1', 'PSME2', 'XAF1', 'IRF7', 'IRF9', 'AIM2', 'ISG20']
    
    #            'AFFX-HUMISGF3A/M97935_3_at','AFFX-HUMISGF3A/M97935_MA_at','AFFX-HUMISGF3A/M97935_MB_at']
    
    #Antigen Recognition genes within first 40 iterations
    SSVM_HLA = ['HLA-DQA1', 'HLA-DQB1', 'HLA-DOB', 'HLA-E', 'HLA-B', 'HLA-DPA1','CD1C',
                'TAP1', 'TAP2','MICA']
    
    #TNF Super Family
    SSVM_TNF = ['TNFRSF10B', 'TNFSF10', 'TNFRSF4', 'TNFAIP1','TNFAIP6', 'TNF', 'TNFAIP3', 'TNFRSF14']
    #'TNFRSF9' at iteration 41
     
    #IL-1 Beta Receptor Family
    SSVM_IL1 = ['IL1B', 'IL1R1','IL1RAP','IL1F5','IL1RL2','IL33']
    
    
    #B Cell Maturation and Activation
    SSVM_BCell = [  'IGHM','IGHD', 'IGHV3-23',
                   'CD200',  'CD24','CD9','CD22', 'CD38', 'CD19', 'CD79A', 'CD86', 'CD40',
                   'CD72' ] 
    
    #Cell cycle related genes
    SSVM_CDC = ['CDK5','CDKAL1', 'CDK5R2','CDCA8','CDC45L','CDKN1C', 'CDC20', 'CDCA3', 'CDKL5' ]
        
    #RNA Helicases DDX series 'RTP4',
    #SSVM_DDX = [ 'DDX60', 'DDX17', 'DDX58']
    
    #Programmed Cell Death
    SSVM_PCD = ['PDCD1LG2' ,'PCDHGA11', 'PDCD4', 'PCDHA3', 'CASP7','CASP5','CASP4','CASP10' ]
    
    #Chemokines
    SSVM_CHEMO = ['CX3CR1','CXCL6','CXCL11','CXCR5','CXCL10','CCL5','CCL11', 'CCR1',
                      'CCR3','CCR6','CCR10','CCRL2','DARC']

    #Cell Adhesion Molecules
    SSVM_CAM = ['ICAM4', 'ICAM3', 'ICAM5', 'MADCAM1']
    
    #Cytokine-Cytokine Interaction
    SSVM_CYT = ['IL16', 'IL17RC', 'IL18RAP', 'IL22', 'IL9', 'SOCS1', 'SOCS6', 'SOCS3']
    
    #Other Immune Response
    SSVM_IR = ['KIR2DL3','FCER2','FCER1G','FCRL2','CD8A','LY6E', 'LY9', 'MARCO', 'TLR5']
    
    #Complement Pathway
    SSVM_CP = ['CR2', 'C2', 'C3AR1', 'C1QA', 'C1QB']
    
        
    PATHWAYS = {1:("Interferon Stimulated Genes", SSVM_ISG),
                2:("Antigen Recognition Genes", SSVM_HLA),
                3:("TNF Super Family", SSVM_TNF),
                4:("IL-1 Beta Receptor Family", SSVM_IL1),
                5:("B Cell Maturation and Activation", SSVM_BCell),
                6:("Cell Cycle Related", SSVM_CDC),
                7:("Programmed Cell Death", SSVM_PCD),
                8:("Chemokines", SSVM_CHEMO),
                9:("Cell Adhesion Molecules", SSVM_CAM),
                10:("Cytokine-Cytokine Receptor Signaling", SSVM_CYT),
                11:("Complement Pathway", SSVM_CP),
                12:("Other Immune Response", SSVM_IR) }
    

if __name__ == '__main__':
    pass