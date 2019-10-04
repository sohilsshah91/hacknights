from __future__ import print_function
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, Formula
pandas2ri.activate()
from rpy2.robjects.packages import importr
dexseq = importr('DEXSeq')
bp = importr('BiocParallel')
'''
Adopted from: https://stackoverflow.com/questions/41821100/running-deseq2-through-rpy2
'''

to_dataframe = robjects.r('function(x) data.frame(x)')

class py_DEXSeq:
    '''
    DEXSeq2 object through rpy2
    input:
    count_matrix: should be a pandas dataframe with each column as count, and a id column for exon id
        example:
        exonID    sampleA    sampleB
        geneA    5    1
        geneB    4    5
        geneC    1    2
    design_matrix: an design matrix in the form of pandas dataframe, see DESeq2 manual, samplenames as rownames
                treatment
    sampleA1        A
    sampleA2        A
    sampleB1        B
    sampleB2        B
    design_formula: see DEXSeq manual, example: "~ sample + exon + exon:treatment""
    feature_column: column name of exon id columns, example "id"
    var_column: will pass to fitExpToVar for DEXSeq exon fold change
    exons: exon id
    genes: gene id for dexseq grouping
    threads: number of threads to use
    '''
    def __init__(self, count_matrix, design_matrix, design_formula,
                feature_column='id', var_column = 'condition',
                exons=None, genes=None, threads=1):
        try:
            assert feature_column in count_matrix.columns, 'Wrong gene id column name'
            assert var_column in design_matrix.columns, 'Wrong var column for DEXSeq'
        except AttributeError:
            sys.exit('Wrong Pandas dataframe?')

        self.dxd = None
        self.dxd_res = None
        self.dexseq_result = None
        self.comparison = None
        self.normalized_count_matrix = None
        self.feature_column = feature_column
        self.exons = exons
        self.genes = genes
        self.gene_id = count_matrix[self.feature_column]
        self.count_matrix = pandas2ri.py2ri(count_matrix.drop(feature_column,axis=1))
        self.design_matrix = pandas2ri.py2ri(design_matrix)
        self.design_formula = Formula(design_formula)
        self.BPPARAM = bp.MulticoreParam(workers=threads)
        self.var_column = var_column


    def run_dexseq(self, **kwargs):

        self.dxd = dexseq.DEXSeqDataSet(countData=self.count_matrix,
                                        sampleData=self.design_matrix,
                                        design=self.design_formula,
                                        featureID = self.exons,
                                        groupID = self.genes)
        print('Constructed DXD object')
        self.dxd = dexseq.estimateSizeFactors_DEXSeqDataSet(self.dxd)
        self.dxd = dexseq.estimateDispersions_DEXSeqDataSet(self.dxd, BPPARAM=self.BPPARAM)
        print('Starting DEXSeq test')
        self.dxd = dexseq.testForDEU(self.dxd, BPPARAM=self.BPPARAM)
        self.dxd = dexseq.estimateExonFoldChanges(self.dxd,
                                                fitExpToVar=self.var_column,
                                                BPPARAM=self.BPPARAM)
        print('Finished DEXSeq fold change')

    def get_dexseq_result(self, **kwargs):

        self.dexseq_result = to_dataframe(dexseq.DEXSeqResults(self.dxd), **kwargs)
        self.dexseq_result = pandas2ri.ri2py(self.dexseq_result) ## back to pandas dataframe
        self.dexseq_result['exons'] = self.exons
        self.dexseq_result['genes'] = self.genes
        self.dexseq_result.drop('genomicData', axis=1, inplace=True)

    def normalized_count(self):
        self.normalized_count_matrix = dexseq.counts(self.dxd,normalized=True)
        return normalized_count_matrix
