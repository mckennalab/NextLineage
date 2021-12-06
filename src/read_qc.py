import json
import argparse
import math
import random

parser = argparse.ArgumentParser(description='Create stats on aligned bam files')
parser.add_argument('--pre_post_stats', help='the pre-post-aligned reads count file ',required=True)
parser.add_argument('--output_filter', help='the output file',required=True)

args = parser.parse_args()

out_edits = open(args.output_filter,"w")
out_edits.write("# parent_name: 'pre_post_filtering_rate'\n")
out_edits.write("# description: 'pre_post_filtering_rate'\n")
out_edits.write("# plot_type: 'table'\n")
out_edits.write("# section_name: 'pre_post_filtering_rate'\n")
out_edits.write("sample\tfiltering_rate\n")

for pre_post in args.pre_post_stats.split(","):
    sample = pre_post[:-18]
    lines_file = open(pre_post)
    pre  = int(lines_file.readline().strip().split()[0])
    post = int(lines_file.readline().strip().split()[0])
    out_edits.write(sample + "\t" + str(post/pre) + "\n")

out_edits.close()



