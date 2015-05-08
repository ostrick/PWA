#!/usr/bin/env python
import sys
# Command line args are in sys.argv[1], sys.argv[2] ...
# sys.argv[0] is the script name itself and can be ignored

import matplotlib.pyplot as plt
from math import *



if len(sys.argv) == 2:
    name = sys.argv[1]
    f = open(name) 

else:
    f = open("Mult-eta15a-fixed-E.dat") 



E = []

E0p = []
E1p = []
M1p = []
M1m = []

E2p = []
E2m = []
M2p = []
M2m = []

E3p = []
E3m = []
M3p = []
M3m = []

E4p = []
E4m = []
M4p = []
M4m = []

E5p = []
E5m = []
M5p = []
M5m = []


E6p = []
E6m = []
M6p = []
M6m = []

E7p = []
E7m = []
M7p = []
M7m = []

E8p = []
E8m = []
M8p = []
M8m = []

E9p = []
E9m = []
M9p = []
M9m = []

k = 0
bereich = [1,3,5,7,9,11,13,15,17]

for line in f:
    words = line.split()
    if ((words[0] == 'W') or (words[0] == '(MeV)')):
        continue
        #    print len(words)
        #    print words[0]

    energy = float(words[0])
    E.append(((938.27*938.27 +  2.*938.27*energy)**0.5))

    k=1
    E0p.append([float(words[k]),float(words[k+1])])        
    E1p.append([float(words[k+2]),float(words[k+3])])        
    M1p.append([float(words[k+4]),float(words[k+5])])        
    M1m.append([float(words[k+6]),float(words[k+7])])        

    k=9
    E2p.append([float(words[k]),float(words[k+1])])        
    E2m.append([float(words[k+2]),float(words[k+3])])        
    M2p.append([float(words[k+4]),float(words[k+5])])        
    M2m.append([float(words[k+6]),float(words[k+7])])        

    k=17
    E3p.append([float(words[k]),float(words[k+1])])        
    E3m.append([float(words[k+2]),float(words[k+3])])        
    M3p.append([float(words[k+4]),float(words[k+5])])        
    M3m.append([float(words[k+6]),float(words[k+7])])        

    k=25
    E4p.append([float(words[k]),float(words[k+1])])        
    E4m.append([float(words[k+2]),float(words[k+3])])        
    M4p.append([float(words[k+4]),float(words[k+5])])        
    M4m.append([float(words[k+6]),float(words[k+7])])        

    k=33
    E5p.append([float(words[k]),float(words[k+1])])        
    E5m.append([float(words[k+2]),float(words[k+3])])        
    M5p.append([float(words[k+4]),float(words[k+5])])        
    M5m.append([float(words[k+6]),float(words[k+7])])        

    k=41
    E6p.append([float(words[k]),float(words[k+1])])        
    E6m.append([float(words[k+2]),float(words[k+3])])        
    M6p.append([float(words[k+4]),float(words[k+5])])        
    M6m.append([float(words[k+6]),float(words[k+7])])        

    k=49
    E7p.append([float(words[k]),float(words[k+1])])        
    E7m.append([float(words[k+2]),float(words[k+3])])        
    M7p.append([float(words[k+4]),float(words[k+5])])        
    M7m.append([float(words[k+6]),float(words[k+7])])        

    k=57
    E8p.append([float(words[k]),float(words[k+1])])        
    E8m.append([float(words[k+2]),float(words[k+3])])        
    M8p.append([float(words[k+4]),float(words[k+5])])        
    M8m.append([float(words[k+6]),float(words[k+7])])        

    k=65
    E9p.append([float(words[k]),float(words[k+1])])        
    E9m.append([float(words[k+2]),float(words[k+3])])        
    M9p.append([float(words[k+4]),float(words[k+5])])        
    M9m.append([float(words[k+6]),float(words[k+7])])        


#a = [ e[0] for e in E0p ]

f.close()

#########################################################
f1 = open('E0p.txt', 'w')
f1.write(" W         E0+(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = E0p[i][0]
    Im = E0p[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('E1p.txt', 'w')
f1.write(" W         E1+(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = E1p[i][0]
    Im = E1p[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('M1p.txt', 'w')
f1.write(" W         M1+(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = M1p[i][0]
    Im = M1p[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('M1m.txt', 'w')
f1.write(" W         M1-(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = M1m[i][0]
    Im = M1m[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()
#########################################################

#########################################################
f1 = open('E2p.txt', 'w')
f1.write(" W         E2+(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = E2p[i][0]
    Im = E2p[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('E2m.txt', 'w')
f1.write(" W         E2-(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = E2m[i][0]
    Im = E2m[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('M2p.txt', 'w')
f1.write(" W         M2+(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = M2p[i][0]
    Im = M2p[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('M2m.txt', 'w')
f1.write(" W         M2-(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = M2m[i][0]
    Im = M2m[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()
#########################################################

#########################################################
f1 = open('E3p.txt', 'w')
f1.write(" W         E3+(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = E3p[i][0]
    Im = E3p[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('E3m.txt', 'w')
f1.write(" W         E3-(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = E3m[i][0]
    Im = E3m[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('M3p.txt', 'w')
f1.write(" W         M3+(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = M3p[i][0]
    Im = M3p[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('M3m.txt', 'w')
f1.write(" W         M3-(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = M3m[i][0]
    Im = M3m[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()
#########################################################

#####################################################
f1 = open('E4p.txt', 'w')
f1.write(" W         E4+(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = E4p[i][0]
    Im = E4p[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('E4m.txt', 'w')
f1.write(" W         E4-(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = E4m[i][0]
    Im = E4m[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('M4p.txt', 'w')
f1.write(" W         M4+(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = M4p[i][0]
    Im = M4p[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('M4m.txt', 'w')
f1.write(" W         M4-(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = M4m[i][0]
    Im = M4m[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()
#########################################################

#########################################################
f1 = open('E5p.txt', 'w')
f1.write(" W         E5+(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = E5p[i][0]
    Im = E5p[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('E5m.txt', 'w')
f1.write(" W         E5-(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = E5m[i][0]
    Im = E5m[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('M5p.txt', 'w')
f1.write(" W         M5+(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = M5p[i][0]
    Im = M5p[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('M5m.txt', 'w')
f1.write(" W         M5-(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = M5m[i][0]
    Im = M5m[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()
#########################################################

#########################################################
f1 = open('E6p.txt', 'w')
f1.write(" W         E6+(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = E6p[i][0]
    Im = E6p[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('E6m.txt', 'w')
f1.write(" W         E6-(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = E6m[i][0]
    Im = E6m[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('M6p.txt', 'w')
f1.write(" W         M6+(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = M6p[i][0]
    Im = M6p[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('M6m.txt', 'w')
f1.write(" W         M6-(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = M6m[i][0]
    Im = M6m[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()
#########################################################

#########################################################
f1 = open('E7p.txt', 'w')
f1.write(" W         E7+(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = E7p[i][0]
    Im = E7p[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('E7m.txt', 'w')
f1.write(" W         E7-(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = E7m[i][0]
    Im = E7m[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('M7p.txt', 'w')
f1.write(" W         M7+(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = M7p[i][0]
    Im = M7p[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('M7m.txt', 'w')
f1.write(" W         M7-(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = M7m[i][0]
    Im = M7m[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()
#########################################################

#########################################################
f1 = open('E8p.txt', 'w')
f1.write(" W         E8+(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = E8p[i][0]
    Im = E8p[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('E8m.txt', 'w')
f1.write(" W         E8-(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = E8m[i][0]
    Im = E8m[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('M8p.txt', 'w')
f1.write(" W         M8+(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = M8p[i][0]
    Im = M8p[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('M8m.txt', 'w')
f1.write(" W         M8-(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = M8m[i][0]
    Im = M8m[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()
#########################################################

#########################################################
f1 = open('E9p.txt', 'w')
f1.write(" W         E9+(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = E9p[i][0]
    Im = E9p[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('E9m.txt', 'w')
f1.write(" W         E9-(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = E9m[i][0]
    Im = E9m[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('M9p.txt', 'w')
f1.write(" W         M9+(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = M9p[i][0]
    Im = M9p[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()

f1 = open('M9m.txt', 'w')
f1.write(" W         M9-(eta_p) \n")
f1.write(" (MeV)     Re      Im \n")
for i in range(len(E)):
    Re = M9m[i][0]
    Im = M9m[i][1]
    f1.write("%4.1f %4.4f %4.4f \n" %(E[i], Re, Im))
f1.close()
#########################################################

