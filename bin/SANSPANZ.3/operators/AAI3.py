import sys
from myoperator import BlockOperator
import DictServer

# Clone a rowwise operator. Class name and file name (operators/classname.py) must be identical.
class AAI3(BlockOperator) :
    """
Output_1 = SANS result plus lineage, genus, kingdom
Output_2 = taxonomic summary of bidirectional best SANS hits per species: weighted AUC, count, genus, kingdom, coverage, average, median, AAI histogram

Filter Output_2 on coverage range, average range, max number of hits to generate Plotly barcharts.
Filter Output_2 on minimum AAI to generate Plotly scatterplot.
Exclude query-species (top .tax species) from Output_1 and select best hit per query to generate Krona taxonomic profile piechart.

"""
    def __init__(self,glob):
        self.glob=glob
        [self.data,self.summary_data]=glob.use_sheets(['data','summary'])
        [self.qpid_col,self.isquery_col,self.pide_col,self.bits_col,self.species_col]=\
                self.data.use_columns(['qpid','isquery','pide','bits','species'])
        # global parameters used: glob.param['PANZ_MAXHITS'], glob.param.speciesindex
        self.summary_data.create_columns(['qpid','csv'])
        # do sums over input stream
        self.counts={} # key = species, counts iindex = pide/100
        self.total={} # key = species, number of hits
        self.wsum={} # key = species, pide-sum
        self.nprot=0
        self.MAXHITS=self.glob.param['PANZ_MAXHITS']
        # count many-to-one matches
        self.pidcount={} # key = spid, count queries
        self.pidmaxpide={} # key = spid, best bidirectional hits
        self.pid2species={} # key = spid, map to species
        # remote access
        self.REMOTE=self.glob.param['CONN_REMOTE']
        self.DICTHOST=self.glob.param['CONN_HOSTNAME']

    def process(self,block):
        seen={}
        self.data.sort_block(self.bits_col,reverse=True)
        if len(block)>0: self.nprot+=1
        nhits=0
        qpid=''
        data=[""]*len(self.glob.speciesindex)
        for row in block:
                qpid=row[self.qpid_col]
                if row[self.isquery_col]=="1": 
                   qpid=row[self.qpid_col]
                   continue # exclude query
                # only count best hit per species
                species=row[self.species_col].rstrip()
                if species in seen: continue
                seen[species]=1
                nhits+=1
                if nhits>self.MAXHITS: break
                if not species in self.glob.speciesindex: continue
                # pide of top-50 species are stored in row-vector
                data[self.glob.speciesindex[species]]=row[self.pide_col]
        if qpid != '': self.summary_data.append_row([qpid,','.join(data)])
