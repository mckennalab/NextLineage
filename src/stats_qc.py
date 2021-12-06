import json
import argparse
import math
import random
import pandas as pd

parser = argparse.ArgumentParser(description='Create stats on aligned bam files')
parser.add_argument('--stats_file', help='The stats file',required=True)
parser.add_argument('--output_edits', help='the output file',required=True)
parser.add_argument('--output_stats', help='the output file',required=True)

args = parser.parse_args()

out_stats = open(args.output_stats,"w")
out_stats.write("# parent_name: 'sample_passing_proportion'\n")
out_stats.write("# description: 'sample_passing_proportion common target reports'\n")
out_stats.write("# plot_type: 'table'\n")
out_stats.write("# section_name: 'Pass-fail counts'\n")
out_stats.write("sample\tPassing_proportion\n")


out_edits = open(args.output_edits,"w")
out_edits.write("# parent_name: 'sample_edit_targets'\n")
out_edits.write("# description: 'sample_edit_targets common target reports'\n")
out_edits.write("# plot_type: 'table'\n")
out_edits.write("# section_name: 'sample_edit_counts'\n")
out_edits.write("sample\ttarget\tevent\tcount\n")

for stats_file in args.stats_file.split(","):
    sample = stats_file[:-9]
    pre_read_ones = pd.read_csv(stats_file,sep="\t")
    pass_val = len(pre_read_ones[pre_read_ones['keep'].str.contains("PASS")])
    fail_val = len(pre_read_ones[pre_read_ones['keep'].str.contains("FAIL")])
    pass_fail = pass_val/(pass_val+fail_val)
    out_stats.write(sample + "\t" + str(pass_fail) + "\n")
    for i in range(1,11):
        top_10 = pre_read_ones[pre_read_ones['keep'] == 'PASS']['target1'].value_counts()[0:10].to_dict()
        for key,value in top_10.items():
            out_edits.write(sample + "\t" + str(i) + "\t" + str(key) + "\t" + str(value) + "\n")
out_edits.close()
out_stats.close()