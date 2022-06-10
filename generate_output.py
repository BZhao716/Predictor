
import os 
import numpy as np 
from Bio.Blast import NCBIXML
import sys

filename = sys.argv[1]
path = sys.argv[2]

blast_dic = os.path.join(path, "programs", "ncbi-blast-2.13.0+", "bin")
blast_db = os.path.join(path, "programs", "ncbi-blast-2.13.0+", "dataset")
disprot_annotation = os.path.join(path, "programs", "ncbi-blast-2.13.0+", "annotation")
compilefile = os.path.join(path, filename)

#run sequence similarity
compileoutput = os.path.join(path, filename+".blast")
formatout = os.path.join(path, filename+".format")
os.chdir(blast_dic)
os.system("./blastp -db "+blast_db+"/DisProt_all_seq.txt -query "+compilefile+" -outfmt '5' -out "+compileoutput)

result = open(compileoutput, "r")
records = NCBIXML.parse(result)
item=next(records)
with open(formatout, "wt") as y:
	for alignment in item.alignments:
		for hsp in alignment.hsps:
			if hsp.expect < 0.000001:
				id = alignment.title.split()[-1]
				align_len = hsp.align_length
				gaps = hsp.gaps
				qstart = hsp.query_start
				qend = hsp.query_end
				sstart = hsp.sbjct_start
				send = hsp.sbjct_end
				qreg = hsp.query
				sreg = hsp.sbjct
				match = hsp.match
				y.write(id+"\t"+str(align_len)+"\t"+str(gaps)+"\t"+str(qstart)+"\t"+str(qend)+"\t"+str(sstart)+"\t"+str(send)+"\n"
					+qreg+"\n"+sreg+"\n"+match+"\n")


#apply similarity to prediction
predout = os.path.join(path, filename+".pred")
disprot_seq = os.path.join(disprot_annotation, "disprot_annotation.txt")
final_output = os.path.join(path, filename+".output")

with open(predout, "rt") as x:
	a = x.readlines()
	query_pred = []
	for i in a:
		query_pred.append(float(i.strip()))

z = open(disprot_seq, "rt")
b = z.read().split(">")[1:]
dic_seq = {}
dic_annot = {}
count = 0
for m in b:
	dic_seq[m.split("\n")[0].strip()] = list(m.split("\n")[1].strip())
	dic_annot[m.split("\n")[0].strip()] = list(map(int,list(m.split("\n")[2].strip())))

q = open(compilefile, "rt")
query = q.readlines()
query_id = query[0].strip()
query_seq = list(query[1].strip())
r = open(formatout, "rt")
rt = r.readlines()
print (len(query_seq), len(query_pred))
with open(final_output, "wt") as u:
	if len(rt)==0:
		u.write(query_id+"\n")
		query_binary = np.where(np.array(query_pred)>=0.3393, 1, 0)
		for uu in range(len(query_seq)):
			u.write(str(uu+1)+"\t"+query_seq[uu]+"\t"+str(query_pred[uu])+"\t"+str(query_binary[uu])+"\n")
	else:
		u.write(query_id+"\n")
		for nn in range(0, len(rt), 4):
			sbjec_id = rt[nn].split("\t")[0]
			query_start = int(rt[nn].split("\t")[3])
			query_end = int(rt[nn].split("\t")[4])
			sbjec_start = int(rt[nn].split("\t")[5])
			sbjec_end = int(rt[nn].split("\t")[6])
			sbjec_annot = list(np.where(np.array(dic_annot[sbjec_id])==3, 1, 0))
			query_reg = rt[nn+1][:-1]
			sbjec_reg = rt[nn+2][:-1]
			match_resi = rt[nn+3][:-1]
			query_seq_reg = query_seq[query_start-1:query_end]
			query_pred_reg = query_pred[query_start-1:query_end]
			sbjec_annot_reg = sbjec_annot[sbjec_start-1:sbjec_end]
			for t1 in range(len(query_reg)):
				if query_reg[t1] == "-":
					query_pred_reg.insert(t1, "-")
					query_seq_reg.insert(t1, "-")
			for t2 in range(len(sbjec_reg)):
				if sbjec_reg[t2] == "-":
					sbjec_annot_reg.insert(t2, "-")
			for t in range(len(sbjec_annot_reg)):
				if match_resi[t] == query_seq_reg[t] or match_resi[t]=="+":
					query_pred_reg[t] = np.average([query_pred_reg[t], sbjec_annot_reg[t]])
			new_query_pred_reg = []
			for p in query_pred_reg:
				if p == "-":
					pass
				else:
					new_query_pred_reg.append(p)
			query_pred[query_start-1:query_end] = new_query_pred_reg
		query_binary = np.where(np.array(query_pred)>=0.3393, 1, 0)
		for uu in range(len(query_seq)):
			u.write(str(uu+1)+"\t"+query_seq[uu]+"\t"+str(query_pred[uu])+"\t"+str(query_binary[uu])+"\n")

os.chdir(path)
os.system("rm -r "+filename+".blast")
os.system("rm -r "+filename+".format")
os.system("rm -r "+filename+".hhblits")
os.system("rm -r "+filename+".npy")
os.system("rm -r "+filename+".pred")
os.system("rm -r "+filename+".pssm")
os.system("rm -r "+filename+".vsl2")
os.system("rm -r asaq.pred")
os.system("rm -r "+filename.split(".")[0]+".ss2")

