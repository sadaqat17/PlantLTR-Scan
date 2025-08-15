import sys
from myoperator import BlockOperator

# Clone a rowwise operator. Class name and file name (operators/classname.py) must be identical. 
class SANStopHtaxid(BlockOperator) :
	"""Sort SANS hit list and return top H hits + link taxids"""
	def __init__(self,glob):
		[self.data,self.summary_data]	= glob.use_sheets(["data","summary"])
		[self.bits_col] 				= self.data.use_columns(["bits"])
		#[self.bits_col]=glob.use_sheet("data").use_columns(["bits"])
		self.H=glob.param['SANS_H']
		# this is a composite operator
		[self.a]						= glob.use_operators(['taxid'])

	def process(self,block):
		# keep top H hits
		block.sort(key = lambda x: float(x[self.bits_col]), reverse=True)
		del(block[self.H:])

		# assign taxids
		for row in block:
			self.a.process(row)
			
