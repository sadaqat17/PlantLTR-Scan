from myoperator import BlockOperator
import re

class output_DE(BlockOperator):
        """
        Select one line per DE-cluster with the best quality description.

        Creates cluster_data column 'desc','genename'
        Inputs: data columns 'clusid','desc','FF','status'
                cluster_data column 'clusid','genename'
        """

        def __init__(self, glob):
                [self.data,self.cluster_data]=glob.use_sheets(["data","cluster"])
                [self.clusid_col1,self.desc_col1,self.qual_col,self.status_col,self.genename_col1]=self.data.use_columns(['clusid','desc','FF','DE_status',"genename"])
                [self.desc_col2,self.genename_col2,self.clusid_col2]=self.cluster_data.use_columns(["desc","genename","clusid"])
                self.MAXHITS=glob.param['PANZ_MAXHITS']
		# clusters created by single linkage, unexpected words may link to wrong DE; use DE-frequency to downweight outliers that get high FF
                self.FFweight=0.1 # DE_selection_score = FF *( 1 + FFweight*nDE)

        def process(self,block):
		# count descriptions in hitlist
                count={}
                for row in block:
                        desc=row[self.desc_col1]
                        if desc not in count: count[desc]=0
                        count[desc]+=1
                # remember RM2, best FF description per cluster weighted by description-count
                desc={}
                bestqual={}
                seen={}
                for row in block:
                        clusid=row[self.clusid_col1]
                        if not clusid in desc:
                                bestqual[clusid]=0.0
                                desc[clusid]=''
                        if row[self.status_col]=="False": continue
                        testdesc=row[self.desc_col1]
                        if testdesc in seen: continue
                        qual=float(row[self.qual_col])*(1+self.FFweight*count[testdesc])
                        if qual > bestqual[clusid]:
                                bestqual[clusid]=qual
                                desc[clusid]=row[self.desc_col1]
                # gene names by majority vote
                gncnt={}
                totcnt=0.0
                maxcnt=0
                maxgn=""
                for row in block:
                        if row[self.status_col]=="False": continue
                        gn=row[self.genename_col1].upper()
                        if gn=="": continue
                        # exclude gene symbols with underscore
                        if re.search(r'\w+_\w+',gn): continue
                        if not gn in gncnt: gncnt[gn]=0
                        gncnt[gn]+=1
                        totcnt+=1.0
                        if totcnt>=self.MAXHITS: break
                for gn in gncnt.keys():
                        if gncnt[gn]>maxcnt:
                                maxcnt=gncnt[gn]
                                maxgn=gn
                if maxcnt/(1.0+totcnt) <= 0.5: maxgn="" # require majority
                # save in cluster_data
                for row in self.cluster_data.block:
                        clusid=row[self.clusid_col2]
                        if not clusid in desc: desc[clusid]=""
                        row[self.desc_col2]=desc[clusid]
                        row[self.genename_col2]=maxgn # copy winner to every cluster


