import os
from Bio import SeqIO

fr1 = open("../cov_abundance_coverage_0.5/all.csv","r").read().splitlines()
fw = open("all_abundance_with_circle.csv","w")
bed_id = open("../circle_-q/all_bed.id","r").read().splitlines()
for line1 in fr1:
	strain = line1.split(",")[0]
	contigid = line1.split(",")[2]
	phageid = line1.split(",")[1]
	c_p_id = contigid + "__" + phageid
	start = int(line1.split(",")[3])
	end = int(line1.split(",")[4])
	seq_path="/beegfs/home/lxd/6_prophage/3.induced_fastq_ont/1_induce_id_phaster/fna_new_nofilt_circlator/" + strain + ".fna"
	seq_records = SeqIO.parse(seq_path,"fasta")
	dict_len = {}
	len_path = "./fna_contig_len/" + strain + "_len.txt"
	len_fw = open(len_path,"w")
	for record in seq_records:
		length=len(record.seq)
		dict_len[record.id] = length
		len_fw.write(record.id + "\t" + str(length) + "\n")
	range2=[]
	fr_path = "/beegfs/home/lxd/6_prophage/3.induced_fastq_ont/3_induce_id_phage_abundance_from_flye/circle_eccfinder/ecc_finder-main/eccfinder/" + strain + "/" + strain + ".csv"
	if os.path.exists(fr_path):
		fr2 = open(fr_path,"r").read().splitlines()
		for line2 in fr2:
			contigid2=line2.split()[0]
			start2=int(line2.split()[1]) + 1
			end2=int(line2.split()[2]) + 1
			if contigid2 == contigid and end2 > start and start2 < end:
				x = range(start,end+1)
				y = range(start2,end2+1)
				xs = set(x)
				inter = xs.intersection(y)
				inter_len = len(inter)
				x_per = inter_len/(end - start + 1)
				y_per = inter_len/(end2 - start2 + 1)
				if x_per > 0.2 and y_per > 0.2:
					start_end = str(start2) + "-" + str(end2)
					range2.append(start_end)
	eccfinder_ranges = ";".join(range2)
	contig_end = dict_len[contigid]
	range3=[]
	if strain in bed_id:
		fr_path = "../circle_-q/bam/" + strain + "_circle.bed"
		fr2 = open(fr_path,"r").read().splitlines()
		for line2 in fr2:
			contigid2=line2.split()[0]
			start2=int(line2.split()[1]) + 1
			end2=int(line2.split()[2]) + 1
			split_reads = int(line2.split()[4])
			if contigid2 == contigid and end2 > start and start2 < end and split_reads >= 2:
				x = range(start,end+1)
				y = range(start2,end2+1)
				xs = set(x)
				inter = xs.intersection(y)
				inter_len = len(inter)
				x_per = inter_len/(end - start + 1)
				y_per = inter_len/(end2 - start2 + 1)
				if x_per > 0.2 and y_per > 0.2:
					start_end = str(start2) + "-" + str(end2)
					range3.append(start_end)
	circlemap_range = ";".join(range3)
	fw.write(line1 + "," + eccfinder_ranges + "," + circlemap_range + "," + str(contig_end) + "\n")
