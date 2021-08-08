from bs4 import BeautifulSoup
import numpy as np


eV = 1.16045E+04 ### 1eV = 1.16045E+04 K
Hz = 2.41799E+14 ### 1eV = 2.41799E+14 Hz
filelist = ['5E-04', '1E-03', '3E-03', '1E-02', '3E-02', '1E-01']

name = filelist[0]
with open(name+'.html', encoding = 'utf-8') as f:
	html = f.read()

soup = BeautifulSoup(html, 'html.parser')
Energy = soup.find(class_="text-muted").get_text()
Energy = Energy.split("\n")

j = 0
l = list()

for i, text in enumerate(Energy):
	if "Energy" in text:
		l.append(1)
		l[j] = i
		j = j+1

l.append(1)
l[j] = l[j-1]+(l[j-1]-l[j-2])
print(l)
Tmp = np.array([float(n)*1.0e+03*eV for n in filelist])
R = np.zeros(len(l)-1)
nu = np.zeros(l[1]-1-l[0])

kappa = np.zeros((l[1]-1-l[0], len(filelist), len(l)-1))

for j, name in enumerate(filelist):
	with open(name+'.html', encoding = 'utf-8') as f:
		html = f.read()
	
	soup = BeautifulSoup(html, 'html.parser')
	Energy = soup.find(class_="text-muted").get_text()
	Energy = Energy.split("\n")

	for i in range(len(l)-1):
		R[i] = float(Energy[l[i]].split()[10])/Tmp[j]**1.5
		array = np.array([[float(n) for n in Energy[k].split()] for k in range(l[i]+1, l[i+1])])
		kappa[:, j, i] = array[:, 2]
		nu = array[:, 0]

R = np.log10(R)
Tmp = np.log10(Tmp)

np.savetxt('opacity_table/freqency.txt', nu.reshape(-1, 1)*1.0e+03*Hz, fmt = '%1.4e')

for i in range(l[1]-1-l[0]):
	name = 'opacity_table/opacity_{:05d}.txt'.format(i)
	R_list = R.tolist()
	R_list = ['{:1.3f}'.format(n) for n in R_list]
	header = ' '.join(map(str, R_list))
	np.savetxt(name, np.append(Tmp.reshape(-1, 1), np.log10(kappa[i, :, :]), axis=1), header = header, fmt = '%1.3f', comments = '')
