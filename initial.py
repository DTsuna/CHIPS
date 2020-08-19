import pylab
import mesa_reader as mr
import matplotlib.pyplot as plt
pylab

path = 'hoge1.txt'

h = mr.MesaData('profile54' + '.data')
size = len(h.zone)

m_sun = 1.9884e+33
r_sun = 6.96e+10

with open(path, mode = 'w') as f:
#        f.write('j =' + str(size) + '\n')
#        f.write('grid m_r(Msun) m_r(g)  dmass(g) radius(cm) density(g/cm^3) pressure(erg/cm^3) temperature(K) h1 he3 he4 c12 n14 o16 ne20 mg24 si28 s32 ar36 ca40 ti44 cr48 cr56 fe52 fe54 fe56 ni56\n')
        for n in range(size, 0, -1):
                temp_zone = size - h.zone[n-1] + 1
                temp_mr = h.mass[n-1]
                temp_mr_cgs = temp_mr * m_sun
                if n == size:
                        mr_inner = 0
                else:
                        mr_inner = h.mass[n] * m_sun
                temp_dmass = temp_mr_cgs - mr_inner
                temp_rasius = pow(10, h.logR[n-1]) * r_sun
                temp_density = pow(10, h.logRho[n-1])
                temp_pressure = pow(10 ,h.logP[n-1])
                temp_temperature = pow(10, h.logT[n-1])
                temp_h1 = h.h1[n-1]
                temp_he3 = h.he3[n-1]
                temp_he4 = h.he4[n-1]
                temp_c12 = h.c12[n-1]
                temp_n14 = h.n14[n-1]
                temp_o16 = h.o16[n-1]
                temp_ne20 = h.ne20[n-1]
                temp_mg24 = h.mg24[n-1]
                temp_si28 = h.si28[n-1]
                temp_s32 = h.s32[n-1]
                temp_ar36 = h.ar36[n-1]
                temp_ca40 = h.ca40[n-1]
                temp_ti44 = h.ti44[n-1]
                temp_cr48 = h.cr48[n-1]
                temp_cr56 = h.cr56[n-1]
                temp_fe52 = h.fe52[n-1]
                temp_fe54 = h.fe54[n-1]
                temp_fe56 = h.fe56[n-1]
                temp_ni56 = h.ni56[n-1]

                f.write(str(temp_zone) + ' ' + str(temp_mr) + ' ' + str(temp_mr_cgs) + ' ' + str(temp_dmass) + ' ' + str(temp_rasius) + ' ' + str(temp_density) + ' ' + str(temp_pressure) + ' ' + str(temp_temperature) + ' ' + str(temp_h1) + ' ' + str(temp_he3) +  ' ' + str(temp_he4) +  ' ' + str(temp_c12) +  ' ' + str(temp_n14) + ' ' + str(temp_o16) +  ' ' + str(temp_ne20) +  ' ' + str(temp_mg24) +  ' ' + str(temp_si28) +  ' ' + str(temp_s32) +  ' ' + str(temp_ar36) +  ' ' + str(temp_ca40) +  ' ' + str(temp_ti44) +  ' ' + str(temp_cr48) +  ' ' + str(temp_cr56) +  ' ' + str(temp_fe52) +  ' ' + str(temp_fe54) +  ' ' + str(temp_fe56) +  ' ' + str(temp_ni56) + '\n')
