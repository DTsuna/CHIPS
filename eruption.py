import convert
import subprocess

file_me = 'profile54.data'
file_hydro = 'InitForHydro.txt'
hydroNumMesh = 10000
logscaleRemesh = False

massCutByHand = False # If true, massCutPoint is used. If false, helium core is cutted automatically.
massCutPoint = 1.3 # unit in Msun

subprocess.call(["rm", "f/inclmn.f"])
subprocess.call(["rm", "f/eruptPara.d"])
subprocess.call(["rm", "InitForHydro.txt"])

convert.convertForHydro(file_me, file_hydro, hydroNumMesh, massCutByHand, massCutPoint, logscaleRemesh)


injectedEnergy = 1.5e47
injectDuration = 1e3
time_CSM = 5e0


ScaledByEnvelopeEnergy = True # If enabled, (-1)*(injectedEnergyRate)*(total energy of the envelope) is deposited instead of injectedEnergy.
injectedEnergyRate = 0.3 # around 0.3 is recommended

continueTransfer = False # If true, radiative transfer scheme is activated even after the eruption.

convert.setSnhydParam(hydroNumMesh,time_CSM,injectedEnergy,injectDuration, ScaledByEnvelopeEnergy, injectedEnergyRate, continueTransfer)


subprocess.call(["mkdir", "snhydOutput"])
subprocess.call(["make", "clean"])
subprocess.call("make")

subprocess.call("./runsnhyd")
