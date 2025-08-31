#!/beegfs/home/lxd/miniconda3/bin/python
import sys
from Bio import SeqIO
import numpy as np

file1 = sys.argv[1] #Input genomecov file
name=file1.split(".")[0]
#Read in the file of the sequencing depth
path1 = "/beegfs/home/lxd/6_prophage/3.induced_fastq_ont/3_induce_id_phage_abundance_from_flye/circle_bam_genomecov/bam/" + file1
fr1 = open(path1,"r")
lines1 = fr1.read().splitlines()
dict_cov = {}
for line1 in lines1:
	id1 = line1.split()[0]
	cov1 = int(line1.split()[2])
	dict_cov.setdefault(id1,[]).append(cov1)
#calculate the length of each contig
seq_path="/beegfs/home/lxd/6_prophage/3.induced_fastq_ont/1_induce_id_phaster/fna_new_nofilt_circlator/" + name + ".fna"
seq_records = SeqIO.parse(seq_path,"fasta")
dict_len = {}
for record in seq_records:
	length=len(record.seq)
	dict_len[record.id] = length
#calculate the depth of each contig except phage
phage_path="/beegfs/home/lxd/6_prophage/3.induced_fastq_ont/1_induce_id_phaster/phaster_result_nofilt_circlator_summary/" + name +"_summary.txt"
fr_phage=open(phage_path,"r")
phages=fr_phage.read().splitlines()
dict_contig_cov = {}
for contig_id in dict_len.keys():
	contig_cov=[]
	start=1
	for phage in phages[34:]:
		phage_contig_id = phage.split()[4].split(",")[0]
		if phage_contig_id == contig_id:
			end=int(phage.split()[4].rsplit(":",1)[1].split("-")[0])
			contig_cov.extend(dict_cov[contig_id][start-1:end-1])
			start=int(phage.split()[4].rsplit(":",1)[1].split("-")[1])+1
	end=dict_len[contig_id]
	contig_cov.extend(dict_cov[contig_id][start-1:end])
	dict_contig_cov[contig_id]=contig_cov
#calculate the depth of each phage
dict_phage_cov={}
for phage in phages[34:]:
	phage_contig_id = phage.split()[4].split(",")[0]
	phage_score = phage.split()[2]
	start=int(phage.split()[4].rsplit(":",1)[1].split("-")[0])
	end=int(phage.split()[4].rsplit(":",1)[1].split("-")[1])
	phage_id = phage.split()[0]
	phage_cov = dict_cov[phage_contig_id][start-1:end]
	contig_cov=dict_contig_cov[phage_contig_id]
	num_cov=float(0)
	for x in phage_cov:
		if x!=0:
			num_cov=num_cov + 1
	coverage=num_cov/float(len(phage_cov))
	phage_cov_mean = np.mean(phage_cov)
	if len(contig_cov) == 0:
		contig_cov_mean = 0
	else:
		contig_cov_mean = np.mean(contig_cov)
	if phage_cov_mean > contig_cov_mean and coverage > 0.5:
		cov_mean = str(phage_cov_mean - contig_cov_mean)
		dict_phage_cov.setdefault(phage_id,[]).append(phage_contig_id)
		dict_phage_cov.setdefault(phage_id,[]).append(str(start))
		dict_phage_cov.setdefault(phage_id,[]).append(str(end))
		dict_phage_cov.setdefault(phage_id,[]).append(cov_mean)
		dict_phage_cov.setdefault(phage_id,[]).append(str(coverage))
		dict_phage_cov.setdefault(phage_id,[]).append(str(phage_score))
path_w = "./csv/" + name + "_phage_abundance.csv"  
fw1=open(path_w,"w")
cov_mean_sum=float(0)
for phage_id in dict_phage_cov.keys():
	cov_mean_sum = cov_mean_sum + float(dict_phage_cov[phage_id][3])
for phage_id in dict_phage_cov.keys():
	abundance = float(dict_phage_cov[phage_id][3])/cov_mean_sum
	fw1.write(name + "," + "phage" + phage_id + "," + dict_phage_cov[phage_id][0] + "," + dict_phage_cov[phage_id][1] + "," + dict_phage_cov[phage_id][2] + "," + dict_phage_cov[phage_id][3] + "," + dict_phage_cov[phage_id][4] + "," + dict_phage_cov[phage_id][5] + "," + str(abundance) + "\n")
