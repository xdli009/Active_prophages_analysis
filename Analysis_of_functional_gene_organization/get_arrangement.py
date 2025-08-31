
#在位于末端的整合酶家族中，只有一个噬菌体G25-86__p14的整合酶有两个，但是都位于前端，所以不影响本脚本的判断
list_integrase = []
dict_integrase = {}
fr_integrase = open("/beegfs/home/lxd/6_prophage/4.prophage_genetics/att_site/integrase_CDS_pos/integrase_id_subfamily.txt","r").read().splitlines()
for line in fr_integrase:
	l0 = line.split("\t")[0]
	family = line.split("\t")[1]
	if family == "subfamily-duplicated":
		list_integrase.append(l0)
	else:
		phageid = l0.split("_CDS")[0]
		pos = int(l0.split("_CDS_")[1])
		if pos < 3:
			dict_integrase[phageid] = "front"
		else:
			dict_integrase[phageid] = "end"

fr1 = open("id_positive_module_class.txt","r").read().splitlines()
list1 = []
dict_annot = {}
dict_num = {}
dict_label = {}
dict_category = {}
for line in fr1[1:]:
	hl = line.split("\t")[1]
	phageid = line.split("\t")[0]
	path_par = "/beegfs/home/lxd/6_prophage/4.prophage_genetics/2.function_module/3.module_cluster/pharokka_p_n/" + phageid + "_cds_final_merged_output.tsv"
	fr_par = open(path_par,"r").read().splitlines()
	pos = dict_integrase[phageid]
	pre = ""
	if pos == "front":
		for line in fr_par:
			gene_name = line.split("\t")[0]
			label = line.split("\t")[21]
			annot = hl + label
			category = line.split("\t")[22]
			if category == "connector" or category == "DNA, RNA and nucleotide metabolism" or category == "head and packaging" or category == "integration and excision" or category == "lysis" or category == "tail" or category == "transcription regulation":
				if gene_name not in list_integrase:
					if pre != "":
						pairs = pre + "---" + annot
						if pairs in dict_annot.keys():
							dict_annot[pairs] = dict_annot[pairs] + 1
						else:
							dict_annot[pairs] = int(1)
					dict_num.setdefault(annot,[]).append(gene_name)
					dict_category[annot] = category
					dict_label[annot] = label
					pre = annot
	else:
		for line in reversed(fr_par):
			gene_name = line.split("\t")[0]
			label = line.split("\t")[21]
			annot = hl + label
			category = line.split("\t")[22]
			if category == "connector" or category == "DNA, RNA and nucleotide metabolism" or category == "head and packaging" or category == "integration and excision" or category == "lysis" or category == "tail" or category == "transcription regulation":
				if gene_name not in list_integrase:
					if pre != "":
						pairs = pre + "---" + annot
						if pairs in dict_annot.keys():
							dict_annot[pairs] = dict_annot[pairs] + 1
						else:
							dict_annot[pairs] = int(1)
					dict_num.setdefault(annot,[]).append(gene_name)
					dict_category[annot] = category
					dict_label[annot] = label
					pre = annot
fw_path = "id_positive_module_input.txt"
fw = open(fw_path,"w")
fw.write("annot1" + "\t" + "annot2" + "\t" + "num" + "\n")
for pairs in dict_annot.keys():
	annot1 = pairs.split("---")[0]
	annot2 = pairs.split("---")[1] 
	fw.write(annot1 + "\t" + annot2 + "\t" + str(dict_annot[pairs]) + "\n")
fw2_path = "id_positive_module_annot_num.txt"
fw2 = open(fw2_path,"w")
fw2.write("annot" + "\t" + "size" + "\n")
for annot in dict_num.keys():
	num = len(dict_num[annot])
	fw2.write(annot + "\t" + str(num) + "\n")
fw3_path = "id_positive_module_category.txt"
fw3 = open(fw3_path,"w")
fw3.write("annot" + "\t" + "type" + "\n")
for annot in dict_category.keys():
	fw3.write(annot + "\t" + dict_category[annot] + "\n")
fw4_path = "id_positive_module_label.txt"
fw4 = open(fw4_path,"w")
fw4.write("annot" + "\t" + "label" + "\n")
for annot in dict_label.keys():
	fw4.write(annot + "\t" + dict_label[annot] + "\n")
