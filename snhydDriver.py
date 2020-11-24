import convert
import subprocess

file_me = 'profile.data'
file_hydro = 'InitForHydro.txt'
hydroNumMesh = 10000

massCutByHand = False # If true, massCutPoint is used. If false, helium core is cutted automatically.
massCutPoint = 1.3 # unit in Msun

subprocess.call(["rm", "f/inclmn.f"])
subprocess.call(["rm", "f/eruptPara.d"])
subprocess.call(["rm", "InitForHydro.txt"])

convert.convertForHydro(file_me, file_hydro, hydroNumMesh, massCutByHand, massCutPoint)


injectedEnergy = 1.5e47
injectDuration = 1e3
time_CSM = 5e0


ScaledByEnvelopeEnergy = True # If enabled, (-1)*(injectedEnergyRate)*(total energy of the envelope) is deposited instead of injectedEnergy.
injectedEnergyRate = 0.3 # around 0.3 is recommended

convert.setSnhydParam(hydroNumMesh,time_CSM,injectedEnergy,injectDuration, ScaledByEnvelopeEnergy, injectedEnergyRate)


subprocess.call(["mkdir", "snhydOutput"])
subprocess.call(["make", "clean"])
subprocess.call("make")

subprocess.call("./runsnhyd")
