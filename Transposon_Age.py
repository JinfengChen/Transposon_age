#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
from Bio.AlignIO.MafIO import MafIndex
from Bio import AlignIO

def usage():
    test="name"
    message='''
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = unit[1]
    return data

def maf_index(infile, chrs):
    index = MafIndex('%s.index' %(infile), infile, chrs)
    return index

#Chr1    RepeatMasker    Transposon      1004    1114    236     +       .       ID=TE000001;Target=(CCCTAA)n 5 102;Class=Simple_repeat;PercDiv=0.0;PercDel=0.0;PercIns=11.7;
#Chr1    RepeatMasker    Transposon      5684    5763    374     +       .       ID=TE000005;Target=STOWAWAY29_OS 1 83;Class=DNA/TcMar-Stowaway;PercDiv=16.2;PercDel=3.8;PercIns=0.0;
def gff_parse(infile):
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd: 
            line = line.rstrip()
            if len(line) > 2 and not line.startswith('#'): 
                unit  = re.split(r'\t',line)
                start = int(unit[3]) 
                end   = int(unit[4])
                chro  = unit[0]
                strand= unit[6]
                temp  = defaultdict(str)
                attrs = re.split(r';', unit[8])
                for attr in attrs:
                    #print attr
                    if not attr == '':
                        #print 'yes'
                        idx, value = re.split(r'\=', attr)
                        temp[idx] = value
                repid   = temp['ID']
                repname = temp['Target']
                repfam  = temp['Class']
                data[chro].append([repid, start, end, repname, repfam])
                #print '%s\t%s\t%s\t%s\t%s\t%s' %(repid, chro, start, end, repname, repfam)
                #print '%s\t%s\t%s' %(chro, start, end)
    return data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-m', '--maf')
    parser.add_argument('-g', '--gff')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.maf) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = 'transposon.age.OPU.list'

    gff = gff_parse(args.gff)
    ofile = open(args.output, 'w')
    for chro in sorted(gff.keys()):
        chro_id = 'MSU7.%s' %(chro)
        maf     = '%s/%s.axt.chain.prenet.net.axt.maf' %(args.maf, chro)
        index   = maf_index(maf, chro_id)
        #repid, start, end, repname, repfam
        print chro
        for repeat in gff[chro]:
            flag    = 0
            starts  = [repeat[1]]
            ends    = [repeat[2]]
            print repeat[0], repeat[1], repeat[2]
            try:
                multi_align = index.get_spliced(starts, ends, strand = '+1')
                #AlignIO.write(multi_align, '%s.%s.align.maf' %(chro, repeat[0]), 'maf')
                #print repeat[0]
                #conserved in alignment
                if len(multi_align) > 1:
                    ref_len = 0
                    tar_len = 0
                    for seq in multi_align:
                        print '>'
                        print seq.id
                        print seq.seq
                        if str(seq.id) == chro_id:
                            ref_len = len(seq.seq) - seq.seq.count('-')
                        else:
                            ##deal with same chromosome or not?
                            temp_len = len(seq.seq) - seq.seq.count('-')
                            tar_len = temp_len if temp_len > tar_len else tar_len
                    #print '%s\t%s' %(ref_len, tar_len)
                    #cover > 80% of element
                    if float(tar_len)/float(ref_len) > 0.8:
                        flag = 1
                print >> ofile, '%s\t%s\t%s\t%s' %(repeat[0], flag, repeat[2]-repeat[1]+1, repeat[4])
            except:
                print 'can not slice alignment'
                print repeat[0], repeat[1], repeat[2]
                continue
    ofile.close()
if __name__ == '__main__':
    main()

