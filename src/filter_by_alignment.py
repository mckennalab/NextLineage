import subprocess
import argparse
import gzip
import itertools 

parser = argparse.ArgumentParser(description='filter reads against the genome')
parser.add_argument('--reference_genome', help='the whole genome were going to align against',required=True)
parser.add_argument('--reference_gestalt_name', help='what we call the GESTALT contig inside this genome',required=True)
parser.add_argument('--read1', help='the first read (10x) containing the actual transcript sequence',required=True)
parser.add_argument('--barcode', help='barcode, or barcode and UMI for v2/v3 chemistry',required=True)
parser.add_argument('--umi', help='the umi read file')
parser.add_argument('--cpus', help='the number of CPUs',type=int)

parser.add_argument('--output_sample_prefix', help='the prefix for output files (plus read.fq.gz, barcode.fq.gz, and umi.fq.gz',required=True)

args = parser.parse_args()

# why python doesn't have this like any normal lang... https://alexwlchan.net/2018/12/iterating-in-fixed-size-chunks/
def chunked_iterable(iterable, size):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, size))
        if not chunk:
            break
        yield chunk

temp_sam = "temp.sam"
command_array = ["bwa","mem","-t", str(args.cpus), "-o",temp_sam,args.reference_genome,args.read1]
print(" ".join(command_array))
subprocess.run(command_array)

aligned_read_names = {}
samfile = open(temp_sam)
for line in samfile:
    if not line.startswith("@"):
        if args.reference_gestalt_name in line: # catching both primary and secondary alignments 
            aligned_read_names["@" + line.split("\t")[0]] = True

reads    = chunked_iterable(gzip.open(args.read1,'rt'),4)
barcodes = chunked_iterable(gzip.open(args.barcode,'rt'),4)

reads_output    = gzip.open(args.output_sample_prefix + ".read.fq.gz", 'wt')
reads_filtered  = gzip.open(args.output_sample_prefix + ".read_filtered.fq.gz", 'wt')
barcodes_output = gzip.open(args.output_sample_prefix + ".barcode.fq.gz", 'wt')


read_count = 0
collected_reads = 0

if args.umi:
    umis     = chunked_iterable(gzip.open(args.umi,'rt'),4)
    umis_output     = gzip.open(args.output_sample_prefix + ".umi.fq.gz", 'wt')
    for read in reads:
        barcode = next(barcodes)
        umi = next(umis)
        read_count += 1

        if read[0].split(" ")[0] in aligned_read_names:
            reads_output.write(read[0].strip() + '\n' + read[1].strip() + '\n+\n' + read[3].strip() + '\n')
            barcodes_output.write(barcode[0].strip() + '\n' + barcode[1].strip() + '\n+\n' + barcode[3].strip() + '\n')
            umis_output.write(umi[0].strip() + '\n' + umi[1].strip() + '\n+\n' + umi[3].strip() + '\n')
            collected_reads += 1
        else:
            reads_filtered.write(read[0].strip() + '\n' + read[1].strip() + '\n+\n' + read[3].strip() + '\n')
else:
    for read in reads:
        barcode = next(barcodes)
        read_count += 1

        if read[0].split(" ")[0] in aligned_read_names:
            reads_output.write(read[0].strip() + '\n' + read[1].strip() + '\n+\n' + read[3].strip() + '\n')
            barcodes_output.write(barcode[0].strip() + '\n' + barcode[1].strip() + '\n+\n' + barcode[3].strip() + '\n')
            collected_reads += 1


reads_output.close()
reads_filtered.close()
barcodes_output.close()
if args.umi:
    umis_output.close()

print(collected_reads/read_count)