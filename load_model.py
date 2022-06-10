
import tensorflow as tf 
from tensorflow import keras
from tensorflow.keras import layers
import numpy as np
from keras.regularizers import l2, l1 
from keras.models import Model
from keras.models import Sequential
from keras.models import load_model
import sys
from base64 import b64encode
import os

filename = sys.argv[1]
path = sys.argv[2]

model = tf.keras.models.load_model('./model/binding_5-2_model_220525.h5')
query_input = filename + ".npy"

def load_data(input_file):
	data = np.load(input_file)
	return data

def separate_data(dataset):
	X1 = np.ndarray.astype(dataset[:,:,:129], dtype="float")
	X22 = np.ndarray.astype(dataset[:,:,129:-5], dtype="float")
	X33 = np.ndarray.astype(dataset[:,:,-5:], dtype="float")
	X2 = []
	X3 = []
	for i in range(len(X1)):
		X2.append(X22[i][2])
		X3.append(X33[i][2])
	X2 = np.array(X2)
	X3 = np.array(X3)
	print (X1.shape)
	print (X2.shape)
	print (X3.shape)
	return X1, X2, X3

data = load_data(query_input)
x1, x2, x3 = separate_data(data)

with tf.device('/cpu:10'):
	result = model.predict([x2, x1, x3])

t = list(map(str,list(result)))
print (len(list(result)))
print (len(t))

compileoutput = os.path.join(path, filename+".pred")
with open(compileoutput, "w") as x:
	for j in t:
		x.write(j[1:-1]+"\n")

x.close()





 

