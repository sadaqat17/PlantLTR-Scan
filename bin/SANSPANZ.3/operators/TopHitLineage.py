import sys
from myoperator import BlockOperator

# Clone a rowwise operator. Class name and file name (operators/classname.py) must be identical. 
class TopHitLineage(BlockOperator) :
    """Extract top bitscore hits and add lineage."""
    def __init__(self,glob):
        [self.qpid_col,self.bits_col,self.pide_col,self.isquery_col,self.species_col,self.lineage_col]=glob.use_sheet("data").use_columns(["qpid","bits","pide","isquery","species","lineage"])
	self.summary_data=glob.use_sheet("summary")
	self.summary_data.use_columns(["#query","bitscore","identity","lineage","species"])
	[self.l]=glob.use_operators(['lineage'])

    def process(self,block):
	if len(block) == 0: return
	bestbits=0.0
	datarow=[]
#	block.sort(key = lambda x: float(x[self.bits_col]), reverse=True)
	for row in block:
		x=float(row[self.bits_col])
		if x > bestbits:
			bestbits=x
			datarow=row
	# copy first row to summary sheet
#	datarow=block[0]
	if len(datarow)==0: return
	if datarow[self.isquery_col] == '1': 
		lineage=''
	else:
		self.l.process(datarow)
		lineage=datarow[self.lineage_col]
	self.summary_data.append_row([datarow[self.qpid_col],datarow[self.bits_col],datarow[self.pide_col],lineage,datarow[self.species_col]])	
