#!/usr/bin/python
#This script is going to smooth the data by taking average of 11 values
#This script is for bending angle
#Han Wen

import sys
import numpy as np

file_name = str(sys.argv[1])
file_inA = open(file_name+"_A.csv", "r")
file_inB = open(file_name+"_B.csv", "r")
file_inC = open(file_name+"_C.csv", "r")
file_inD = open(file_name+"_D.csv", "r")

file_out = open(file_name+"smooth10.csv", "w")

frame = []
p_inA = []
a_inA = []
p_inB = []
a_inB = []
p_inC = []
a_inC = []
p_inD = []
a_inD = []


lines = file_inA.readlines()
for line in lines:
    words = line.split(',')
    frame.append(float(words[0]))
    p_inA.append(float(words[1]))
    a_inA.append(float(words[2]))

file_inA.close()


lines = file_inB.readlines()
for line in lines:
    words = line.split(',')
    p_inB.append(float(words[1]))
    a_inB.append(float(words[2]))

file_inB.close()


lines = file_inC.readlines()
for line in lines:
    words = line.split(',')
    p_inC.append(float(words[1]))
    a_inC.append(float(words[2]))

file_inC.close()


lines = file_inD.readlines()
for line in lines:
    words = line.split(',')
    p_inD.append(float(words[1]))
    a_inD.append(float(words[2]))

file_inD.close()




for i in range(len(frame)):
    if i<5 or i > len(frame)-6:
	s = str(frame[i]) + "," + str(p_inA[i])+"," +str(a_inA[i]) +"," + str(p_inB[i])+"," +str(a_inB[i])+ "," +str(p_inC[i])+"," +str(a_inC[i])+"," +str(p_inD[i])+"," +str(a_inD[i]) +"\n"
	file_out.write(s)
    else:


	total_tempA = 0#p_inA[i-2] + p_inA[i-1] +p_inA[i] +p_inA[i+1] +p_inA[i+2]
        total_temaA = 0#a_inA[i-2] + a_inA[i-1] +a_inA[i] +a_inA[i+1] +a_inA[i+2]
        total_tempB = 0#p_inB[i-2] + p_inB[i-1] +p_inB[i] +p_inB[i+1] +p_inB[i+2]
        total_temaB = 0#a_inB[i-2] + a_inB[i-1] +a_inB[i] +a_inB[i+1] +a_inB[i+2]
        total_tempC = 0#p_inC[i-2] + p_inC[i-1] +p_inC[i] +p_inC[i+1] +p_inC[i+2]
        total_temaC = 0#a_inC[i-2] + a_inC[i-1] +a_inC[i] +a_inC[i+1] +a_inC[i+2]
        total_tempD = 0#p_inD[i-2] + p_inD[i-1] +p_inD[i] +p_inD[i+1] +p_inD[i+2]
        total_temaD = 0#a_inD[i-2] + a_inD[i-1] +a_inD[i] +a_inD[i+1] +a_inD[i+2]


	for j in range(i-5,i+6):
	    total_tempA += p_inA[j]
            total_temaA += a_inA[j]
            total_tempB += p_inB[j]
            total_temaB += a_inB[j]
            total_tempC += p_inC[j]
            total_temaC += a_inC[j]
            total_tempD += p_inD[j]
            total_temaD += a_inD[j]
   

        ave_pA = total_tempA/10.0
        ave_aA = total_temaA/10.0
        ave_pB = total_tempB/10.0
        ave_aB = total_temaB/10.0
        ave_pC = total_tempC/10.0
        ave_aC = total_temaC/10.0
        ave_pD = total_tempD/10.0
        ave_aD = total_temaD/10.0

        s = str(frame[i]) + "," + str(ave_pA)+"," +str(ave_aA) +"," + str(ave_pB)+"," +str(ave_aB) +"," + str(ave_pC)+"," +str(ave_aC) +"," + str(ave_pD)+"," +str(ave_aD) + "\n"
        file_out.write(s)
file_out.close()

