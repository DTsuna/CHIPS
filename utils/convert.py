# include package
import math
import sys
import gc
import numpy as np
import mesa_reader as mr

def f(CommonR, FirstTerm, SumMass, NumOfTerms):
        return np.power(CommonR,NumOfTerms) - SumMass*CommonR/FirstTerm + SumMass/FirstTerm - 1

def fx(CommonR, FirstTerm, SumMass, NumOfTerms, DeltaDiff):
        return ( f(CommonR+DeltaDiff, FirstTerm, SumMass, NumOfTerms) - f(CommonR, FirstTerm, SumMass, NumOfTerms) ) / DeltaDiff


def getCommonRatio(CommonR_init, FirstTerm, SumMass, NumOfTerms):
        print('get common ratio for logscale remeshing')
        CommonR = CommonR_init
        Error = 1.0
        MaxError = 1e-8
        DeltaDiff = 1e-10
        i = 0
        print('Trial='+str(i)+'  Common Ratio='+str(CommonR)+'  Error='+str(Error))

        while Error > MaxError:
                pref = f(CommonR, FirstTerm, SumMass, NumOfTerms)
                CommonR = CommonR - (f(CommonR, FirstTerm, SumMass, NumOfTerms)/fx(CommonR, FirstTerm, SumMass, NumOfTerms, DeltaDiff))
                Error = np.absolute( f(CommonR, FirstTerm, SumMass, NumOfTerms) - pref )
                i = i + 1
                print('Trial='+str(i)+'  Common Ratio='+str(CommonR)+'  Error='+str(Error))
                if i > 1000:
                        print('logscale remeshing does not converge')
                        sys.exit(1)
        print('========== Consistency check ==========')
        print('Required total mass ='+str(SumMass))
        result = ( np.power(CommonR,NumOfTerms)*FirstTerm - FirstTerm )/( CommonR - 1)
        print('Computed total mass ='+str(result))
        if np.absolute(result - SumMass) > 2e33 * 1e-6:
                print('Total mass after remeshing is inconsistent with before remeshing!')
                sys.exit(1)
        return CommonR

def convertForHydro(inputFile, outputFile, hydroNumMesh, massCutByHand, massCutPoint, logscaleRemesh):

        #################################################################
        ##################### Input Parameters ##########################
        #################################################################
        # input file from MESA
        h = mr.MesaData(inputFile)
        originalSize = len(h.zone)

        # output file for hydro
        path = outputFile
        size = hydroNumMesh # number of mesh
        if massCutByHand == False:
                massCut = h.he_core_mass + 0.2 # In solar mass unit
        if massCutByHand == True:
                massCut = massCutPoint

        print('massCut =' + str(massCut))
        if logscaleRemesh == True:
                print('use logscale remesh')
        if logscaleRemesh == False:
                print('use linear remesh')

        lowerLimX = 1.0e-40 # Lower limit on the value of composition X
        elemNum = 19 # number of element

        # physical constant
        MSUN = 1.9884e+33
        RSUN = 6.96e+10
        G = 6.6743e-8
        #################################################################
        #################################################################
        #################################################################
        print('Remeshing process is proceeding ...')

        # read data from mesa
        originalMr = np.zeros(originalSize)
        originalMrCgs = np.zeros(originalSize)
        originalDmass = np.zeros(originalSize)
        originalRadius = np.zeros(originalSize)
        originalDensity = np.zeros(originalSize)
        originalPressure = np.zeros(originalSize)
        originalTemperature = np.zeros(originalSize)
        originalX = np.zeros((elemNum, originalSize))
        missing_elem = []

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
                try:
                    originalX[10][i] = h.ar36[originalSize - i - 1]
                except AttributeError:
                    originalX[10][i] = lowerLimX
                    missing_elem.append('ar36') if 'ar36' not in missing_elem else None
                try:
                    originalX[11][i] = h.ca40[originalSize - i - 1]
                except AttributeError:
                    originalX[11][i] = lowerLimX
                    missing_elem.append('ca40') if 'ca40' not in missing_elem else None
                try:
                    originalX[12][i] = h.ti44[originalSize - i - 1]
                except AttributeError:
                    originalX[12][i] = lowerLimX
                    missing_elem.append('ti44') if 'ti44' not in missing_elem else None
                try:                    
                    originalX[13][i] = h.cr48[originalSize - i - 1]
                except AttributeError:
                    originalX[13][i] = lowerLimX
                    missing_elem.append('cr48') if 'cr48' not in missing_elem else None
                try:
                    originalX[14][i] = h.cr56[originalSize - i - 1]
                except AttributeError:
                    originalX[14][i] = lowerLimX
                    missing_elem.append('cr56') if 'cr56' not in missing_elem else None
                try:
                    originalX[15][i] = h.fe52[originalSize - i - 1]
                except AttributeError:
                    originalX[15][i] = lowerLimX
                    missing_elem.append('fe52') if 'fe52' not in missing_elem else None
                try:
                    originalX[16][i] = h.fe54[originalSize - i - 1]
                except AttributeError:
                    originalX[16][i] = lowerLimX
                    missing_elem.append('fe54') if 'fe54' not in missing_elem else None
                try:
                    originalX[17][i] = h.fe56[originalSize - i - 1]
                except AttributeError:
                    originalX[17][i] = lowerLimX
                    missing_elem.append('fe56') if 'fe56' not in missing_elem else None
                try:
                    originalX[18][i] = h.ni56[originalSize - i - 1]
                except AttributeError:
                    originalX[18][i] = lowerLimX
                    missing_elem.append('ni56') if 'ni56' not in missing_elem else None

        if len(missing_elem) > 0:
                print('Warning: Attribute %s is not found and thus ignored...' % ', '.join(missing_elem)) 

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

        if logscaleRemesh == False:
                fiducialDmass = (cuttedMrCgs[originalSize - cellCut - 1] - cuttedMrCgs[0])/(size - 1)
                dmass[0] = cuttedMrCgs[0];
                for i in range(1, size):
                        dmass[i] = ((-1.6/(size-2))*(i-1) + 1.8) * fiducialDmass
                for i in range(0, size):
                        mrCgs[i] = 0
                        for j in range(0, i + 1):
                                mrCgs[i] = mrCgs[i] + dmass[j]

        if logscaleRemesh == True:
                dmass[0] = cuttedMrCgs[0];
                CommonR_init = 1.01
                FirstTerm = 4.0*3.14*cuttedRadius[originalSize - cellCut - 1]*cuttedRadius[originalSize - cellCut - 1]*0.1
                SumMass = cuttedMrCgs[originalSize - cellCut - 1] - cuttedMrCgs[0]
                NumOfTerms = size-1
                CommonR = getCommonRatio(CommonR_init, FirstTerm, SumMass, NumOfTerms)
                print('Frist term ='+str(FirstTerm)+'  Common ratio ='+str(CommonR))
                dmass[size - 1] = FirstTerm
                for i in range(2, size):
                        dmass[size - i] = dmass[size - i + 1]*CommonR
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
        for i in range(1,size):
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


def setSnhydParam(hydroNumMesh,timeToCC,injectedEnergy,injectDuration, ScaledByEnvelopeEnergy, injectedEnergyRate, continueTransfer):
        if ScaledByEnvelopeEnergy == True:
                flag = 1
        else:
                flag = 0
        if continueTransfer == True:
                flag2 = 1
        else:
                flag2 = 0
        with open('src/eruption/f/inclmn.f', mode = 'w') as f:
                f.write('      integer mn, nelem\n')
                f.write('      parameter ( mn = '+str(hydroNumMesh+ 10)+', nelem = 19 )\n')
        with open('src/eruption/f/eruptPara.d', mode = 'w') as f2:
                f2.write('TimeToCC InjectedEnergy InjectDuration ScaledByEnvelopeEnergy injectedEnergyRate\n')
                f2.write(str(timeToCC*86400*365.25) + ' ' +  str(injectedEnergy) + ' ' +  str(injectDuration) + ' ' + str(flag) + ' ' + str(injectedEnergyRate) + ' ' + str(flag2) + '\n')





