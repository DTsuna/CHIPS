import numpy as np
import os

def fopen_op(filename):
	with open(filename) as f:
		R = f.readline().strip('\n').split(' ')
	
	a = np.loadtxt(filename, skiprows = 1)
	temp = [r[0] for r in a]
	kappa = np.delete(a, 0, 1)
	return R, np.array(temp), kappa


def gen_op_frq(Y, dir_name):
	array_Y = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
	array_Y_dir_name = ['Y{:1.1f}'.format(n) for n in array_Y]
	size = len(array_Y)
	for i in range(size):
		if Y >= array_Y[i] and Y < array_Y[i+1]:
			break
	
	freq_file = os.path.dirname(__file__)+'/frequency.txt'
	freq = np.loadtxt(freq_file)
	freq_num = len(freq)
	
	for j in range(freq_num):
		filename = os.path.dirname(__file__)+'/'+array_Y_dir_name[i]+'/opacity_table/opacity_{:05d}.txt'.format(j)
		filename1 = os.path.dirname(__file__)+'/'+array_Y_dir_name[i+1]+'/opacity_table/opacity_{:05d}.txt'.format(j)
		R, temp, kappa = fopen_op(filename)
		R, temp, kappa1 = fopen_op(filename1)
		kappa_outp = ((array_Y[i+1]-Y)*kappa+(Y-array_Y[i])*kappa1)/(array_Y[i+1]-array_Y[i])

		temp = temp.reshape(-1, 1)
		X = np.append(temp, kappa_outp, axis = 1)

		header = ' '.join(R)
		opacity_file = dir_name+'/opacity_{:05d}.txt'.format(j)
		np.savetxt(opacity_file, X, header = header, comments = '', fmt = '%1.4f')
