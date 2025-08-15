from __future__ import absolute_import
from __future__ import print_function
#!/usr/bin/python

import Parameters, Runner
import sys
#sys.path.append('/home/luholm/AAI-profiler')

if __name__ == '__main__' :
    example="""
This is a customized script that performs AAI-profiling on protein fasta sequences.

Usage:
        python %s -i proteinfastafile [-R]

Outputs:
        tax.tab = onesided AAI profiles per species
        tax2.tab = bidirectional AAI profiles per species
        kt.txt = onesided krona input excluding top species
        kt_self.txt = onesided krona input including top species
        matrix.txt = sequence identities of hits to top 50 species

Shell script for AAI-profiler:
        python runsanspanz.py -m SANS -i inputfastafile -o /data/tmp/HID/sans.tab
        python AAI-profiler-tab.py -i /data/tmp/HID/sans.tab -d ./
        perl tax2scatter.pl "onesided AAI-profiler results" onesided ./ 0.5 20 0 1 0 1 < ./tax.tab > ~/public_html/AAI.html
        perl tax2scatter.pl "bidirectional AAI-profiler results" bidirectional ./ 0.5 20 0 1 0 1 < ./tax2.tab > ~/public_html/AAI2.html
        /usr/local/bin/ktImportText -u http://krona.sourceforge.net -o ~/public_html/krona.html ./kt.txt
        /usr/local/bin/ktImportText -u http://krona.sourceforge.net -o ~/public_html/krona_self.html ./kt_self.txt
        perl heatmap.pl > ~/public_html/matrix.html

        """ %sys.argv[0]
    if len(sys.argv) < 3 :
        print(example, file=sys.stderr)
        sys.exit(1)
    try :
        glob=Parameters.WorkSpace()
        USEREMOTE=glob.param['CONN_REMOTE']
        MAXHITS=2000 # include all hits
        # use remote access to SANSparallel and DictServer
        glob.param['CONN_REMOTE'] = USEREMOTE
        print('# Remote access is',USEREMOTE, file=sys.stderr)
        # output files
        sanstabfile='sans.tab'
        taxtabfile='tax.tab'
        taxtabfile2='tax2.tab'
        kttxtfile='kt.txt'
        ktselftxtfile='kt_self.txt'
        matrixtabfile='matrix.txt'
        liveData='liveData'
        # step 0: run SANS
        print('# Running SANS over the Internet. Progress monitored in ',liveData, file=sys.stderr)
        z=Runner.Runner(glob,operator_name='SANS',CHUNK=100,liveData=liveData)
        z.lazyRunner(glob.param['input_FILE'],[sanstabfile],input_format="FASTA")
        # step 1: collect AAI profiles per species
        glob.param['PANZ_MAXHITS']=MAXHITS
        glob.param['input_QUERYSPECIES']='species'
        glob.kt_format=False
        print('# Running AAI2. MAXHITS is ',glob.param['PANZ_MAXHITS'], file=sys.stderr)
        z=Runner.Runner(glob,operator_name='AAI',liveData=None)
        z.lazyRunner(sanstabfile, [None,taxtabfile], input_format='tab', block_column_name='qpid' )
        # step 2: get best hit from sans.tab including top species
        # re-initialize workspace to get clean sheet headers
        glob=Parameters.WorkSpace()
        glob.param['CONN_REMOTE'] = USEREMOTE
        glob.param['PANZ_MAXHITS']=1
        glob.param['input_QUERYSPECIES']='species'
        glob.kt_format=True
        print('# Extracting best hit. All species included. MAXHITS is ',glob.param['PANZ_MAXHITS'], file=sys.stderr)
        z=Runner.Runner(glob,operator_name='AAI',CHUNK=1000,liveData=None)
        z.lazyRunner(sanstabfile, [None,ktselftxtfile], input_format='tab', block_column_name='qpid' )
        # step 3: get best hit from sans.tab excluding top species
        glob=Parameters.WorkSpace()
        glob.param['CONN_REMOTE'] = USEREMOTE
        glob.param['PANZ_MAXHITS']=1
        glob.kt_format=True
        fh=open('tax.tab','r')
        species_col=fh.readline().rstrip().split("\t").index('species')
        glob.param['input_QUERYSPECIES']=fh.readline().split("\t")[species_col].rstrip()
        fh.close()
        print('# Extracting best hit. %s excluded. MAXHITS is ' %glob.param['input_QUERYSPECIES'],glob.param['PANZ_MAXHITS'], file=sys.stderr)
        z=Runner.Runner(glob,operator_name='AAI',CHUNK=1000,liveData=None)
        z.lazyRunner(sanstabfile, [None,kttxtfile], input_format='tab', block_column_name='qpid' )
        # * bidirectional AAI *
        glob=Parameters.WorkSpace()
        glob.param['PANZ_MAXHITS']=MAXHITS
        glob.param['input_QUERYSPECIES']='species'
        glob.kt_format=False
        print('# Running AAI2. MAXHITS is ',glob.param['PANZ_MAXHITS'], file=sys.stderr)
        z=Runner.Runner(glob,operator_name='AAI2',liveData=None)
        z.lazyRunner(sanstabfile, [None,taxtabfile2], input_format='tab', block_column_name='qpid' )
        # step 4: output matrix of ordered query proteins vs. pide of top-50 species
        glob=Parameters.WorkSpace()
        glob.param['CONN_REMOTE'] = USEREMOTE
        glob.param['PANZ_MAXHITS']=MAXHITS
        glob.param['input_QUERYSPECIES']='species'
        glob.kt_format=False
        # save top-50 species from kt.txt in dict
        H=50
        glob.speciesindex={}
        fh=open(taxtabfile,'r')
        i=0
        for line in fh.readlines()[0:H]:
            x=line.split('\t')
            species=x[6]
            if species=='species': continue
            glob.speciesindex[species]=i
            i+=1
        fh.close()
        print('# Creating pide matrix of top-%i species' %H, file=sys.stderr)
        z=Runner.Runner(glob,operator_name='AAI3',liveData=None) # uses local dict glob.speciesindex
        z.lazyRunner(sanstabfile, [None,matrixtabfile], input_format='tab', block_column_name='qpid')

    except KeyboardInterrupt :
        print("CLIENT> killed by user...\n", file=sys.stderr)
        sys.exit(1)

