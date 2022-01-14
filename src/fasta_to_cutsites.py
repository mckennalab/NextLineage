import argparse
from os.path import exists
import re 
import collections

def read_fasta(fasta_file):
  sequence = ""
  fasta = open(fasta_file)
  header = fasta.readline()
  assert(header[0] == ">")
  for line in fasta:
    sequence += line.strip()
  return(sequence)

def read_sites_file(sites_file):
  crispr_sites = []
  fasta = open(sites_file)
  header = fasta.readline().strip()
  assert(header == "sites")
  for line in fasta:
    crispr_sites.append(line.strip().upper())
    print(line.strip().upper())
  return(crispr_sites)

def comp(base):
  if base.upper() == "A":
    return("T")
  elif base.upper() == "T":
    return("A")
  elif base.upper() == "C":
    return("G")
  elif base.upper() == "G":
    return("C")
  else:
    return("N")

def rev_comp(string):
  rev = [comp(x) for x in string]
  rev.reverse()
  return("".join(rev))

class CRISPRInfo:
  def __init__(self,length,cutsite):
    self.length = length
    self.cutsite = cutsite

def cas_to_target_length_and_cutsite_location(cas_type):
  if cas_type.upper() == "CAS12A":
    return(CRISPRInfo(24,20))
  elif cas_type.upper() == "CAS9":
    return(CRISPRInfo(23,17))
  else:
    raise NameError("Unknown CRISPR type:" + cas_type)

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='Create lineage pipeline files from a reference, primers, and known cutsites.')
  parser.add_argument('--reference', help='The reference file. All auxiliary files will be created with the this path as as a basename',required=True)
  parser.add_argument('--sites', help='a file, with a header "stats", containing a single target site per line afterwards. No order is assumed, and multiple hits can be found within the reference for any one target',required=True)
  parser.add_argument('--forward_primer', help='the forward primer, which can be used in the lineage pipeline to ensure the read originates from this target site',required=True)
  parser.add_argument('--reverse_primer', help='the reverse primer, which can be used in the lineage pipeline to ensure the read originates from this target site',required=True)
  parser.add_argument('--CRISPR_type', help='The CRISPR enzyme class in use. For instance, cas9, cas12a, cas9_abe, etc',required=True)
  
  args = parser.parse_args()
  
  crispr_length_cutsite = cas_to_target_length_and_cutsite_location(args.CRISPR_type)

  primers_file  = args.reference + ".primers"
  cutsites_file = args.reference + ".cutSites"
  
  reference_sequence = read_fasta(args.reference).upper()
  sites              = read_sites_file(args.sites)

  if exists(primers_file) or exists(cutsites_file):
    raise NameError("Unable to overwrite existing primers or cutsites file")

  # output the primers file once we've checked that they exist in the reference
  if (not args.forward_primer in reference_sequence) or (not args.reverse_primer in reference_sequence):
    raise NameError("Unable to find the forward or reverse primer in the reference sequence")
  output_primers  = open(primers_file,"w")
  output_primers.write(args.forward_primer + "\n" + args.reverse_primer + "\n")
  output_primers.close()

  # now look for cutsites in the reference
  output_cutsites = open(cutsites_file, "w")
  output_cutsites.write("sites\tposition\tcutPos\n")

  found_sites = []
  for site in sites:
    rev_comp_str = rev_comp(site)

    fwd_result = [_.start() for _ in re.finditer(site, reference_sequence)] 
    rev_result = [_.start() for _ in re.finditer(rev_comp_str, reference_sequence)] 

    for i in fwd_result:
      found_sites.append([i,site,True])

    for i in rev_result:
      found_sites.append([i,site,False])

  sorted_found_sites = sorted(found_sites, key=lambda hit: hit[0])   # sort by position
  
  for hit in sorted_found_sites:
    cutsite = hit[0] + crispr_length_cutsite.cutsite
    if not hit[2]:
      cutsite = hit[0] + (crispr_length_cutsite.length - crispr_length_cutsite.cutsite)

    output_cutsites.write(str(hit[1]) + "\t" + str(hit[0] + 1) + "\t" + str(cutsite + 1) + "\n")

  output_cutsites.close()


 
