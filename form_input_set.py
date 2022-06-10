
import sys
import numpy as np
import pandas as pd
import os
import tempfile
import json

filename = sys.argv[1]
path = sys.argv[2]
input_numpy = filename + ".npy"

#the AA_index used to normalize ASA raw data
max_value = {"A":129,
             "R":274,
             "N":195,
             "D":193,
             "C":167,
             "E":223,
             "Q":225,
             "G":104,
             "H":224,
             "I":197,
             "L":201,
             "K":236,
             "M":224,
             "F":240,
             "P":159,
             "S":155,
             "T":172,
             "W":285,
             "Y":263,
             "V":174}

#the aa_index used to calculate PSSM and obtain the position specific conservation
blosum62_bg = { 'A' : 7.4,
                'R' : 5.2,
                'N' : 4.5,
                'D' : 5.3,
                'C' : 2.5,
                'Q' : 3.4,
                'E' : 5.4,
                'G' : 7.4,
                'H' : 2.6,
                'I' : 6.8,
                'L' : 9.9,
                'K' : 5.8,
                'M' : 2.5,
                'F' : 4.7,
                'P' : 3.9,
                'S' : 5.7,
                'T' : 5.1,
                'W' : 1.3,
                'Y' : 3.2,
                'V' : 7.3  }

bbg_AAs = list(blosum62_bg.keys())
tot = sum(blosum62_bg.values())
blosum62_bg = np.array([ blosum62_bg[k] / tot for k in blosum62_bg ])

#residue index based on DisProt annotation
disorder_index =  {"A":0.0965, "C":-0.5033, "D":0.1880, "E":0.3061, "F":-0.3207, "G":0.1310, "H":-0.1280, "I":-0.3476, "K":0.1817, "L":-0.2850, 
			       "M":-0.1593, "N":-0.0145, "P":0.3119, "Q":0.1447, "R":-0.0220, "S":0.2679, "T":0.0256, "V":-0.2370, "W":-0.4515, "Y":-0.3092}
binding_index = {"A":0.0972, "C":-0.5041, "D":0.2059, "E":0.3440, "F":-0.2785, "G":0.1267, "H":-0.1028, "I":-0.3698, "K":0.2115, "L":-0.2843, 
			     "M":-0.2095, "N":-0.0509, "P":0.3452, "Q":0.1135, "R":-0.0246, "S":0.2210, "T":0.0112, "V":-0.1982, "W":-0.4489, "Y":-0.3793}
function_index = {"A":0.1103, "C":-0.4937, "D":0.1866, "E":0.3212, "F":-0.2977, "G":0.1624, "H":-0.1565, "I":-0.3565, "K":0.1968, "L":-0.3055, 
			      "M":-0.1705, "N":-0.0059, "P":0.3169, "Q":0.1016, "R":-0.0320, "S":0.2719, "T":0.0108, "V":-0.2196, "W":-0.4825, "Y":-0.3487}

def format_psipred(path, id):
	filename = id.split(".")[0] + ".ss2"
	compilename = os.path.join(path, filename)
	b = np.loadtxt(compilename, skiprows=2, dtype="str")
	H = np.array(list(map(float, b[:,4])))
	E = np.array(list(map(float, b[:,5])))
	C = np.array(list(map(float, b[:,3])))
	scoreH = H/(H+E+C)
	scoreE = E/(H+E+C)
	scoreC = C/(H+E+C)
	return scoreH, scoreE, scoreC

def format_asaquick(path):
	filename = "asaq.pred"
	compilename = os.path.join(path, filename)
	data = np.loadtxt(compilename, dtype="str")
	seq = data[:,1]
	scores = data[:,2]
	value = []
	for i in range(len(scores)):
		if seq[i] == "A":
			v = float(scores[i])/max_value["A"]
			if v > 1:
				value.append("1")
			else:
				value.append(str("{:.3f}".format(v)))
		elif seq[i] == "R":
			v = float(scores[i])/max_value["R"]
			if v > 1:
				value.append("1")
			else:
				value.append(str("{:.3f}".format(v)))
		elif seq[i] == "N":
			v = float(scores[i])/max_value["N"]
			if v > 1:
				value.append("1")
			else:
				value.append(str("{:.3f}".format(v)))
		elif seq[i] == "D":
			v = float(scores[i])/max_value["D"]
			if v > 1:
				value.append("1")
			else:
				value.append(str("{:.3f}".format(v)))
		elif seq[i] == "C":
			v = float(scores[i])/max_value["C"]
			if v > 1:
				value.append("1")
			else:
				value.append(str("{:.3f}".format(v)))
		elif seq[i] == "E":
			v = float(scores[i])/max_value["E"]
			if v > 1:
				value.append("1")
			else:
				value.append(str("{:.3f}".format(v)))
		elif seq[i] == "Q":
			v = float(scores[i])/max_value["Q"]
			if v > 1:
				value.append("1")
			else:
				value.append(str("{:.3f}".format(v)))
		elif seq[i] == "G":
			v = float(scores[i])/max_value["G"]
			if v > 1:
				value.append("1")
			else:
				value.append(str("{:.3f}".format(v)))
		elif seq[i] == "H":
			v = float(scores[i])/max_value["H"]
			if v > 1:
				value.append("1")
			else:
				value.append(str("{:.3f}".format(v)))
		elif seq[i] == "I":
			v = float(scores[i])/max_value["I"]
			if v > 1:
				value.append("1")
			else:
				value.append(str("{:.3f}".format(v)))
		elif seq[i] == "L":
			v = float(scores[i])/max_value["L"]
			if v > 1:
				value.append("1")
			else:
				value.append(str("{:.3f}".format(v)))
		elif seq[i] == "K":
			v = float(scores[i])/max_value["K"]
			if v > 1:
				value.append("1")
			else:
				value.append(str("{:.3f}".format(v)))
		elif seq[i] == "M":
			v = float(scores[i])/max_value["M"]
			if v > 1:
				value.append("1")
			else:
				value.append(str("{:.3f}".format(v)))
		elif seq[i] == "F":
			v = float(scores[i])/max_value["F"]
			if v > 1:
				value.append("1")
			else:
				value.append(str("{:.3f}".format(v)))
		elif seq[i] == "P":
			v = float(scores[i])/max_value["P"]
			if v > 1:
				value.append("1")
			else:
				value.append(str("{:.3f}".format(v)))
		elif seq[i] == "S":
			v = float(scores[i])/max_value["S"]
			if v > 1:
				value.append("1")
			else:
				value.append(str("{:.3f}".format(v)))
		elif seq[i] == "T":
			v = float(scores[i])/max_value["T"]
			if v > 1:
				value.append("1")
			else:
				value.append(str("{:.3f}".format(v)))
		elif seq[i] == "W":
			v = float(scores[i])/max_value["W"]
			if v > 1:
				value.append("1")
			else:
				value.append(str("{:.3f}".format(v)))
		elif seq[i] == "Y":
			v = float(scores[i])/max_value["Y"]
			if v > 1:
				value.append("1")
			else:
				value.append(str("{:.3f}".format(v)))
		elif seq[i] == "V":
			v = float(scores[i])/max_value["V"]
			if v > 1:
				value.append("1")
			else:
				value.append(str("{:.3f}".format(v)))
		elif seq[i]== "X" or seq[i]=="U" or seq[i]=="B" or seq[i]=="O" or seq[i]=="Z":
			value.append("0")
	return value

#function to calculate PSSM
def calc_pssm_freq_df(pssm_df):
    df = pssm_df[bbg_AAs]
    freq = np.exp(df) * blosum62_bg
    #freq_sum = np.sum(freq,axis=1)
    return np.array(freq)/np.sum(freq, axis=1)[:,None]

def calc_relative_entropy(pssm_df):
    df = calc_pssm_freq_df(pssm_df)
    cons = np.sum(df*np.log(df/blosum62_bg), axis=1)
    return cons

def format_mmseq(path, id):
	filename = id + ".pssm"
	compilename = os.path.join(path, filename)
	pssm_df = pd.read_csv(compilename, sep="\t", skiprows=1)
	res = calc_relative_entropy(pssm_df)
	score = []
	for t in res:
		score.append(str("{:.2f}".format(t)))
	return pssm_df.values, score

def format_vsl2b(path, id):
	filename = id + ".vsl2"
	compilename = os.path.join(path, filename)
	x = open(compilename, "rt")
	a = x.read().split("Prediction Scores")[1]
	lst = a.split("\n")[4:-2]
	score = []
	for i in lst:
		score.append(float(i.split()[2].strip()))
	return score

def format_hhblits(path, id):
	filename = id + ".hhblits"
	data = np.loadtxt(filename, delimiter="\t", dtype="str")
	return data

def format_index(seq, dic):
	lst = []
	for i in seq:
		if i == "X":
			lst.append(dic["A"])
		elif i == "B":
			lst.append(dic["N"])
		elif i == "Z":
			lst.append(dic["Q"])
		elif i == "U":
			lst.append(dic["M"])
		elif i == "O":
			lst.append(dic["P"])
		else: 	
			lst.append(dic[i])
	return lst

x = open(filename, "rt")
a = x.readlines()
seq = a[1].strip()

#form residue level input
psipredh, psiprede, psipredc = format_psipred(path, filename)
asaquick = format_asaquick(path)
pssm, mmseq = format_mmseq(path, filename)
vsl2 = format_vsl2b(path, filename)
hhblits = format_hhblits(path, filename)
dis = format_index(seq, disorder_index)
bind = format_index(seq, binding_index)
func = format_index(seq, function_index)

#form region level input
ws = 8
new_psiprede = [0]*ws + list(psiprede) + [0]*ws
new_psipredc = [0]*ws + list(psipredc) + [0]*ws
new_asaquick = [0]*ws + list(map(float,asaquick)) + [0]*ws
new_vsl2 = [0]*ws + vsl2 + [0]*ws
new_mmseq = [0]*ws + list(map(float,mmseq)) + [0]*ws
new_dis = [0]*ws + dis + [0]*ws
new_bind = [0]*ws + bind + [0]*ws
new_func = [0]*ws + func + [0]*ws
psiprede_reg = []
psipredc_reg = []
asaquick_reg = []
vsl2_reg = []
mmseq_reg = []
dis_reg = []
bind_reg = []
func_reg = []

for i in range(len(new_psiprede)-16):
	psiprede_reg.append(-np.average(new_psiprede[i:i+17]))
	psipredc_reg.append(np.average(new_psipredc[i:i+17]))
	asaquick_reg.append(np.average(new_asaquick[i:i+17]))
	vsl2_reg.append(np.average(new_vsl2[i:i+17]))
	mmseq_reg.append(-max(new_mmseq[i:i+17]))
	dis_reg.append(np.average(new_dis[i:i+17]))
	bind_reg.append(np.average(new_bind[i:i+17]))
	func_reg.append(np.average(new_func[i:i+17]))


#form protein level input 
psipredc_ind = np.where(psipredc>=np.percentile(psipredc, 90))[0]
psipredc_lst = []
for m in psipredc_ind:
	psipredc_lst.append(psipredc[m])
psipredc_prot = [np.average(psipredc_lst)]*len(psipredc)

vsl2_ind = np.where(vsl2>np.percentile(vsl2, 90))[0]
vsl2_lst = []
for n in vsl2_ind:
	vsl2_lst.append(vsl2[n])
vsl2_prot = [np.average(vsl2_lst)]*len(seq)

dis_prot = [np.average(dis)]*len(seq)
bind_prot = [np.average(bind)]*len(seq)
func_prot = [np.average(func)]*len(seq)

#y = open(input_file, "wt")
ll = []
for j in range(len(seq)):
	ll.append(str(psipredh[j])+"\t"+str(psiprede[j])+"\t"+str(psipredc[j])+"\t"+str(asaquick[j])
		+"\t"+str(vsl2[j])+"\t"+str(mmseq[j])+"\t"+("\t".join(list(map(str,pssm[j][2:]))))+"\t"+str(dis[j])
		+"\t"+str(bind[j])+"\t"+str(func[j])+"\t"+"\t".join(hhblits[j])+"\t"+str(psiprede_reg[j])+"\t"+str(psipredc_reg[j])
		+"\t"+str(asaquick_reg[j])+"\t"+str(vsl2_reg[j])+"\t"+str(mmseq_reg[j])+"\t"+str(dis_reg[j])
		+"\t"+str(bind_reg[j])+"\t"+str(func_reg[j])+"\t"+str(psipredc_prot[j])+"\t"+str(vsl2_prot[j])
		+"\t"+str(dis_prot[j])+"\t"+str(bind_prot[j])+"\t"+str(func_prot[j])+"\n")
	#y.write(str(psipredh[j])+"\t"+str(psiprede[j])+"\t"+str(psipredc[j])+"\t"+str(asaquick[j])
	#	+"\t"+str(vsl2[j])+"\t"+str(mmseq[j])+"\t"+("\t".join(list(map(str,pssm[j][2:]))))+"\t"+str(dis[j])
	#	+"\t"+str(bind[j])+"\t"+str(func[j])+"\t"+"\t".join(hhblits[j])+"\t"+str(psiprede_reg[j])+"\t"+str(psipredc_reg[j])
	#	+"\t"+str(asaquick_reg[j])+"\t"+str(vsl2_reg[j])+"\t"+str(mmseq_reg[j])+"\t"+str(dis_reg[j])
	#	+"\t"+str(bind_reg[j])+"\t"+str(func_reg[j])+"\t"+str(psipredc_prot[j])+"\t"+str(vsl2_prot[j])
	#	+"\t"+str(dis_prot[j])+"\t"+str(bind_prot[j])+"\t"+str(func_prot[j])+"\n")


dd = []
for m in range(len(ll)):
	dd.append(np.array(ll[m].strip().split("\t")))
ws_reg = list(np.ndarray.astype(np.zeros((2, len(dd[0]))), dtype="str"))
new_data = ws_reg + list(dd) + ws_reg
lst = []
for jj in range(len(new_data)-4):
	lst.append(new_data[jj:jj+5])
print (len(lst))
np.save(input_numpy, lst)









