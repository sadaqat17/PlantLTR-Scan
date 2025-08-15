import sys
from myoperator import BlockOperator

# Clone a rowwise operator. Class name and file name (operators/classname.py) must be identical. 
class SANStopH(BlockOperator) :
    """Sort SANS hit list and return top H hits"""
    def __init__(self,glob):
        [self.bits_col]=glob.use_sheet("data").use_columns(["bits"])
        self.H=glob.param['SANS_H']

    def process(self,block):
        block.sort(key = lambda x: float(x[self.bits_col]), reverse=True)
        del(block[self.H:])
