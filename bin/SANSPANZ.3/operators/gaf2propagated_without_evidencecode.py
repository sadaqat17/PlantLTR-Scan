# cut -f 2,5,7 goa_uniprot_noIEA.gaf

from myoperator import RowOperator
from Read_and_Print import read_dict_data
import sys,re,math

class gaf2propagated_without_evidencecode(RowOperator):
    def __init__(self,glob):
        # set parameters
        self.glob=glob
        # propagate assigned GOids to parents
        obotabfile='obo.tab'
        self.directparents=read_dict_data(obotabfile,'goid','directParentSet',header=True,UseColNames=True) # goid,parentset
        # write blockwise summmary
        [self.data,self.godict_data]=glob.use_sheets(["data","godict"])
        # column indices in GAF file
        [self.goid_col1]=self.data.use_columns(["goid"])
        self.godict_data.use_columns(["qoid","propagated"])

    def process(self,row):
        goidlist=[]
        x={}
        GOweight={} # initialize for each block
        goid=row[self.goid_col1]
        goid=re.sub(r'GO:','',goid)
        # direct assignment + parents
        x[goid]=1
        # propagate to all parents, max evidence weight
        stack=x.keys()
        parents={}
        while len(stack)>0:
            goid=stack.pop()
            parents[goid]=1
            if not goid in self.directparents: continue
            tmp=self.directparents[goid]
            if tmp=='': continue # root
            for p in tmp.split(','):
                    if not p in parents: stack.append(p)
        propagated=parents.keys()
        propagated.sort()
        # write to godict spreadsheet
        self.godict_data.append_row([goid,','.join(propagated)])

