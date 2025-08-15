# cut -f 2,5,7 goa_uniprot_noIEA.gaf

from myoperator import BlockOperator
from Read_and_Print import read_dict_data
from PannzerFunctions import Propagate
import sys,re,math

class GOpropagate(BlockOperator):
    """
Input: qpid, goid, soore (FIXED ORDER)
Output: propagated scores

    """
    def __init__(self,glob):
        # set parameters
        self.glob=glob
        self.glob.use_online_dictionaries(['GOIDELIC'])
        # write blockwise summmary
        [self.data,self.summary_data]=glob.use_sheets(["data","summary"])
        self.summary_data.use_columns(['qpid','goid','score'])
        self.data.use_columns(['qpid','goid','score','parents'])
        # column indices in GAF file
        self.qpid_col=0
        self.goid_col=1
        self.score_col=2
        self.parent_col=3

    def process(self,block):
        if len(block)<1: return
        qpid=block[0][self.qpid_col]
        score={}
        for row in block: 
                goid=row[self.goid_col]
                row[self.parent_col]=' '.join(self.glob.GOparents[goid])
                s=float(row[self.score_col])
                score[goid]=s
		# propagate
                if not goid in self.glob.GOparents: continue # root
                for p in self.glob.GOparents[goid]:
                        if p=='': continue
                        if p=='n.d.': continue
                        if not p in score: score[p]=0.0
                        if s>score[p]: score[p]=s
        # output propagated score
        for goid in score:
                datarow=[qpid,goid,str(score[goid])]
                self.summary_data.append_row(datarow)

