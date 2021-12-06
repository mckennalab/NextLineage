import gzip
import argparse
import math
import numpy as np
import argparse
import itertools



parser = argparse.ArgumentParser(description='Join sequencing reads into a single UMIed sequence')
parser.add_argument('--read1', help='the first read (10x) containing the actual transcript sequence',required=True)
parser.add_argument('--barcode', help='the barcode read file (contains the UMI in v2/v3 chemistry)',required=True)
parser.add_argument('--umi', help='the umi read file')
parser.add_argument('--outputfastq', help='the output fastq file of barcode-umi-read',required=True)
parser.add_argument('--trimbases', help='optionally trim off bases from the front of the first read (for instance, hex priming bases)',default=0,type=int)

args = parser.parse_args()

# why python doesn't have this like any normal lang... https://alexwlchan.net/2018/12/iterating-in-fixed-size-chunks/
def chunked_iterable(iterable, size):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, size))
        if not chunk:
            break
        yield chunk

reads = chunked_iterable(gzip.open(args.read1,'rt'),4)
barcodes = chunked_iterable(gzip.open(args.barcode,'rt'),4)
output =  gzip.open(args.outputfastq, 'wt')

if args.umi:
    umis = chunked_iterable(gzip.open(args.umi,'rt'),4)
    for read in reads:
        barcode = next(barcodes)
        umi = next(umis)
        output.write(read[0].strip() + '\n' + barcode[1].strip() + umi[1].strip() + read[1][args.trimbases:len(read[1])].strip() + '\n+\n' + barcode[3].strip() + umi[3].strip() + read[3][args.trimbases:len(read[3])].strip() + '\n')
else:
    for read in reads:
        barcode = next(barcodes)
        output.write(read[0].strip() + '\n' + barcode[1].strip() + read[1][args.trimbases:len(read[1])].strip() + '\n+\n' + barcode[3].strip() + read[3][args.trimbases:len(read[3])].strip() + '\n')

output.close()

