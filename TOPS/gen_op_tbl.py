import numpy as np

size = 8

def fopen_op(filename):
	with open(filename) as f:
		R = f.readline().strip('\n').split(' ')
	
	a = np.loadtxt(filename, skiprows = 1)
	temp = [r[0] for r in a]
	kappa = np.delete(a, 0, 1)
	rowlen = len(temp)
	collen = len(kappa[0])
	return R, np.array(temp), kappa, rowlen, collen


def gen_op_tbl(Y):
	array_Y = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
	array_Y_dir_name = ['Y{:1.1f}'.format(n) for n in array_Y]
	for i in range(size):
		if Y >= array_Y[i] and Y < array_Y[i+1]:
			break
	
	filename = './'+array_Y_dir_name[i]+'/'+array_Y_dir_name[i]+'.txt'
	filename1 = './'+array_Y_dir_name[i+1]+'/'+array_Y_dir_name[i+1]+'.txt'
	R, temp, kappa, rowlen, collen = fopen_op(filename)
	R, temp, kappa1, rowlen, collen = fopen_op(filename1)
	kappa_outp = ((array_Y[i+1]-Y)*kappa+(Y-array_Y[i])*kappa1)/(array_Y[i+1]-array_Y[i])

	temp = temp.reshape(-1, 1)
#	print(np.append(temp, kappa_outp, axis = 1))
	X = np.append(temp, kappa_outp, axis = 1)

	header = ' '.join(R)
	np.savetxt('opacity.txt', X, header = header, comments = '', fmt = '%1.4f')

gen_op_tbl(0.45)
