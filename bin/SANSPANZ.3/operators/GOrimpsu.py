from myoperator import RowOperator
import sys
from Read_and_Print import read_dict_data

class GOrimpsu(RowOperator):
        """
        Add GO class and count vectors to data.

        Creates data columns 'GOclass','GOclass_count'
        Inputs: data columns 'spid','status'
        """
        def __init__(self, glob):
                sys.stderr.write("# Init GOrimpsu %s / %s\n" %(glob.goslimfile,glob.param['input_GOSLIM']))
                self.glob=glob
                [self.go_col,self.gocount_col,self.spid_col,self.status_col]=glob.use_sheet("data").use_columns(['GOclass','GOclass_count','spid','status'])
                glob.use_online_dictionaries(["GOIDELIC","GODICT"])
                self.gofilter = glob.param['input_GOSLIM']
                if self.gofilter:
                       self.golist=read_dict_data(self.gofilter,0,0)
                       sys.stderr.write("#loaded %i GO classes from %s\n" %(len(self.golist),self.gofilter))	

        def process(self,row,verbose=False):
                # no GO terms imported for hits removed by Filter
                #sys.stderr.write("%s\n" %('\t'.join(row)))
                go=""
                cnt=""
                if row[self.status_col] != "False":
                  spid=row[self.spid_col]
                  try:
                    (db,accession,pid)=spid.split("|") # tr|E0SNF6|E0SNF6_DICD3
                    if accession in self.glob.GOdict:
                        tmp=self.glob.GOdict[accession].split(',')
                        #sys.stderr.write("got %s from %s: %s || %s \n" %(accession,spid,accession in self.glob.GOdict, ' '.join(tmp)))
                        for goid in tmp:
                             #sys.stderr.write('goid %s %s %s go %s cnt %s\n' %(goid,goid in self.golist,self.glob.GOcounts[goid],go,cnt))
                             if (not self.gofilter) or (goid in self.golist):
                                 go+=goid+' '
                                 cnt+=self.glob.GOcounts[goid]+" "
                                 #sys.stderr.write("# GOSLIM added goid %s\n" %goid)
                             #else:
                             #    sys.stderr.write("# GOSLIM rejected goid %s\n" %goid)
                  except:
                          if verbose: sys.stderr.write("# Warning: no GOcounts for spid = %s\n" %(spid))
                #row[self.go_col]=" ".join(go).rstrip() # data is string
                row[self.go_col]=go.rstrip()
                row[self.gocount_col]=cnt.rstrip() # data is string

