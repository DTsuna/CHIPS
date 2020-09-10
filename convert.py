# include package
import pylab
import math
import sys
import gc
import numpy as np
import mesa_reader as mr
import matplotlib.pyplot as plt


#################################################################
##################### Input Parameters ##########################
#################################################################
# input file from MESA
h = mr.MesaData('profile54' + '.data')
originalSize = len(h.zone)

# output file for hydro
path = 'output.txt'
size = 10000 # number of mesh
massCut = 4.0 # In solar mass unit
lowerLimX = 1.0e-40 # Lower limit on the value of composition X
elemNum = 19 # number of element

# physical constant
MSUN = 1.9884e+33
RSUN = 6.96e+10
G = 6.6743e-8
#################################################################
#################################################################
#################################################################


# read data from mesa
originalMr = np.zeros(originalSize)
originalMrCgs = np.zeros(originalSize)
originalDmass = np.zeros(originalSize)
originalRadius = np.zeros(originalSize)
originalDensity = np.zeros(originalSize)
originalPressure = np.zeros(originalSize)
originalTemperature = np.zeros(originalSize)
originalX = np.zeros((elemNum, originalSize))

for i in range(0, originalSize):
        originalMr[i] = h.mass[originalSize - i - 1]
        originalMrCgs[i] = originalMr[i] * MSUN
        if i == 0:
                mrInner = 0
        else:
                mrInner = originalMrCgs[i - 1]
        originalDmass[i] = originalMrCgs[i] - mrInner
        originalRadius[i] = pow(10, h.logR[originalSize - i - 1]) * RSUN
        originalDensity[i] = pow(10, h.logRho[originalSize - i - 1])
        originalPressure[i] = pow(10 ,h.logP[originalSize - i - 1])
        originalTemperature[i] = pow(10, h.logT[originalSize - i - 1])                
        originalX[0][i] = h.h1[originalSize - i - 1]
        originalX[1][i] = h.he3[originalSize - i - 1]
        originalX[2][i] = h.he4[originalSize - i - 1]
        originalX[3][i] = h.c12[originalSize - i - 1]
        originalX[4][i] = h.n14[originalSize - i - 1]
        originalX[5][i] = h.o16[originalSize - i - 1]
        originalX[6][i] = h.ne20[originalSize - i - 1]
        originalX[7][i] = h.mg24[originalSize - i - 1]
        originalX[8][i] = h.si28[originalSize - i - 1]
        originalX[9][i] = h.s32[originalSize - i - 1]
        originalX[10][i] = h.ar36[originalSize - i - 1]
        originalX[11][i] = h.ca40[originalSize - i - 1]
        originalX[12][i] = h.ti44[originalSize - i - 1]
        originalX[13][i] = h.cr48[originalSize - i - 1]
        originalX[14][i] = h.cr56[originalSize - i - 1]
        originalX[15][i] = h.fe52[originalSize - i - 1]
        originalX[16][i] = h.fe54[originalSize - i - 1]
        originalX[17][i] = h.fe56[originalSize - i - 1]
        originalX[18][i] = h.ni56[originalSize - i - 1]

"""########################### debug part ###########################
with open(path, mode = 'w') as f:
#        f.write('j =' + str(originalSize) + '\n')
#        f.write('zone m_r(Msun) m_r(g)  dmass(g) radius(cm) density(g/cm^3) pressure(erg/cm^3) temperature(K) h1 he3 he4 c12 n14 o16 ne20 mg24 si28 s32 ar36 ca40 ti44 cr48 cr56 fe52 fe54 fe56 ni56\n')
        for i in range(0, originalSize):
                f.write(str(i + 1) + ' ' + str(originalMr[i]) + ' ' + str(originalMrCgs[i]) + ' ' + str(originalDmass[i]) + ' ' + str(originalRadius[i]) + ' ' + str(originalDensity[i]) + ' ' + str(originalPressure[i]) + ' ' + str(originalTemperature[i]))
                for j in xrange(0,elemNum):
                        f.write(' ' + str(originalX[j][i]))
                f.write('\n')
"""##################################################################


# extract envelope and set lower limit on x
cellCut = originalSize
for i in range(originalSize - 1, -1, -1):
        if originalMr[i] > massCut:
                cellCut = i
if cellCut == originalSize:
        print('EEEOR: Invalid massCut')
        sys.exit(1)

cuttedMrCgs = np.zeros(originalSize - cellCut)
cuttedMr = np.zeros(originalSize - cellCut)
cuttedDmass = np.zeros(originalSize - cellCut)
cuttedRadius = np.zeros(originalSize - cellCut)
cuttedDensity = np.zeros(originalSize - cellCut)
cuttedPressure = np.zeros(originalSize - cellCut)
cuttedTemperature = np.zeros(originalSize - cellCut)
cuttedX = np.zeros((elemNum, originalSize - cellCut))

for i in range(0, originalSize - cellCut):
        cuttedMrCgs[i] = originalMrCgs[i + cellCut]
        cuttedMr[i] = originalMr[i + cellCut]
        cuttedDmass[i] = originalDmass[i + cellCut]
        cuttedRadius[i] = originalRadius[i + cellCut]
        cuttedDensity[i] = originalDensity[i + cellCut]
        cuttedPressure[i] = originalPressure[i + cellCut]
        cuttedTemperature[i] = originalTemperature[i + cellCut]
        for j in range(0, elemNum):
                cuttedX[j, i] = max(lowerLimX, originalX[j, i + cellCut])

del originalMr
del originalMrCgs
del originalDmass
del originalRadius
del originalDensity
del originalPressure
del originalTemperature
del originalX
gc.collect()

"""########################### debug part ###########################
with open(path, mode = 'w') as f:
#        f.write('j =' + str(originalSize) + '\n')
#        f.write('zone m_r(Msun) m_r(g)  dmass(g) radius(cm) density(g/cm^3) pressure(erg/cm^3) temperature(K) h1 he3 he4 c12 n14 o16 ne20 mg24 si28 s32 ar36 ca40 ti44 cr48 cr56 fe52 fe54 fe56 ni56\n')
        for i in range(0, originalSize - cellCut):
                f.write(str(i + 1 + cellCut) + ' ' + str(cuttedMr[i]) + ' ' + str(cuttedMrCgs[i]) + ' ' + str(cuttedDmass[i]) + ' ' + str(cuttedRadius[i]) + ' ' + str(cuttedDensity[i]) + ' ' + str(cuttedPressure[i]) + ' ' + str(cuttedTemperature[i]))
                for j in range(0,elemNum):
                        f.write(' ' + str(cuttedX[j][i]))
                f.write('\n')
"""##################################################################


# remesh profile for hydro
mrCgs = np.zeros(size)
dmass = np.zeros(size)
radius = np.zeros(size)
density = np.zeros(size)
pressure = np.zeros(size)
temperature = np.zeros(size)
x = np.zeros((elemNum, size))

fiducialDmass = (cuttedMrCgs[originalSize - cellCut - 1] - cuttedMrCgs[0])/(size - 1)
dmass[0] = cuttedMrCgs[0];
for i in range(1, size):
        dmass[i] = ((-1.6/(size-2))*(i-1) + 1.8) * fiducialDmass
for i in range(0, size):
        mrCgs[i] = 0
        for j in range(0, i + 1):
                mrCgs[i] = mrCgs[i] + dmass[j]

#print(str(mrCgs[size - 1]))
#print(str(cuttedMrCgs[originalSize - cellCut - 1]))

radius[0] = cuttedRadius[0]
radius[size - 1] = cuttedRadius[originalSize - cellCut -1]
pressure[0] = cuttedPressure[0]
pressure[size - 1] = cuttedPressure[originalSize - cellCut - 1]
temperature[0] = cuttedTemperature[0]
temperature[size - 1] = cuttedTemperature[originalSize - cellCut - 1]
for j in range(0, elemNum):
        x[j][0] = cuttedX[j][0]
        x[j][size - 1] = cuttedX[j][originalSize - cellCut - 1]

logX = np.zeros(elemNum)
for i in range(1, size - 1):
        j = 0
        while mrCgs[i] > cuttedMrCgs[j + 1]:
                j = j + 1

        logRadius = math.log10(cuttedRadius[j]) + (math.log10(cuttedRadius[j+1]) - math.log10(cuttedRadius[j]))*(mrCgs[i] - cuttedMrCgs[j])/(cuttedMrCgs[j+1] - cuttedMrCgs[j])
        logPressure = math.log10(cuttedPressure[j]) + (math.log10(cuttedPressure[j+1]) - math.log10(cuttedPressure[j]))*(mrCgs[i] - cuttedMrCgs[j])/(cuttedMrCgs[j+1] - cuttedMrCgs[j])
        logTemperature = math.log10(cuttedTemperature[j]) + (math.log10(cuttedTemperature[j+1]) - math.log10(cuttedTemperature[j]))*(mrCgs[i] - cuttedMrCgs[j])/(cuttedMrCgs[j+1] - cuttedMrCgs[j])
        for k in range(0, elemNum):
                logX[k] = math.log10(cuttedX[k][j]) + (math.log10(cuttedX[k][j+1]) - math.log10(cuttedX[k][j]))*(mrCgs[i] - cuttedMrCgs[j])/(cuttedMrCgs[j+1] - cuttedMrCgs[j])

        radius[i] = math.pow(10, logRadius)
        pressure[i] = math.pow(10, logPressure)
        temperature[i] = math.pow(10, logTemperature)
        for k in range(0, elemNum):
                x[k][i] = math.pow(10, logX[k])

density[0] = cuttedDensity[0]
for i in range(1, size):
        density[i] = (3.0*(mrCgs[i] - mrCgs[i-1]))/(4.0*math.pi*(math.pow(radius[i],3)-math.pow(radius[i-1],3)))
                
del cuttedMr
del cuttedMrCgs
del cuttedDmass
del cuttedRadius
del cuttedDensity
del cuttedPressure
del cuttedTemperature
del cuttedX
gc.collect()


# shift pressure by half a cell because 
copiedPressure = np.zeros(size)
for i in range(0,size):
        copiedPressure[i] = pressure[i]
for i in xrange(1,size):
        pressure[i] = math.sqrt(copiedPressure[i]*copiedPressure[i-1])


# output
with open(path, mode = 'w') as f:
        f.write(str(size) + '\n')
        f.write('j m_r(Msun) m_r(g) dmass(g) radius(cm) density(g/cm^3) pressure(erg/cm^3) temperature(K) h1 he3 he4 c12 n14 o16 ne20 mg24 si28 s32 ar36 ca40 ti44 cr48 cr56 fe52 fe54 fe56 ni56\n')
        for i in range(0, size):
                f.write(str(i + 1) + ' ' + str(mrCgs[i]/MSUN) + ' ' + str(mrCgs[i]) + ' ' + str(dmass[i]) + ' ' + str(radius[i]) + ' ' + str(density[i]) + ' ' + str(pressure[i]) + ' ' + str(temperature[i]))
                for j in range(0,elemNum):
                        f.write(' ' + str(x[j][i]))
                f.write('\n')











