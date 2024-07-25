import numpy as np  
import scipy as sp
import matplotlib.pyplot as plt
import random as r
import math as mt


#Setting:
L=20 #size of the lattice
beta=0.2 #inverse temperature

Lattice=[]
for i in range(L**2):
    if  r.random()<0.5:
        Lattice.append(1)   
    else:
        Lattice.append(-1)


#ALTERNATIVE: if i wanted to use matrix instead of list, cool for visualizing
Lattice2=[]
for i in range(L**2):
    Lattice2.append(r.randrange(-1, 2, 2))
Lattice2 = np.matrix(Lattice2).reshape((L, L))
#print(Lattice2)



#Mapping to the matrix (labelled from 0)
def line(n):
    return mt.floor(n/L)

def column(n):
    return n%L


#Hamiltonian
def LocalHamiltonian(i,Latticetemp):
    if  i==0:
        return  -Latticetemp[i]*(Latticetemp[i+1]+Latticetemp[i+L]+Latticetemp[i+L-1]+Latticetemp[L*(L-1)] )
    if  i==(L-1):
        return  -Latticetemp[i]*(Latticetemp[i-1]+Latticetemp[0]+Latticetemp[i+L]+Latticetemp[L*L-1])
    if  i==L*(L-1):
        return  -Latticetemp[i]*(Latticetemp[i+1]+Latticetemp[L*L-1]+Latticetemp[0]+Latticetemp[L*(L-2)])
    if  i==L*L-1:
        return  -Latticetemp[i]*(Latticetemp[L*L-2]+Latticetemp[L*(L-1)]+Latticetemp[L-1]+ Latticetemp[L*(L-1)-1])
    if  i%L==0 and i!=0 and i!=L*(L-1):
        return  -Latticetemp[i]*(Latticetemp[i+1]+Latticetemp[i+L-1]+Latticetemp[i-L]+Latticetemp[i+L] )
    if  mt.floor(i/L)==0 and i!=(L-1) and i!=0:
        return  -Latticetemp[i]*(Latticetemp[i-1]+Latticetemp[i+1]+Latticetemp[i+L]+Latticetemp[L*(L-1)+i%L])
    if  mt.floor(i/L)==(L-1) and i!=L*(L-1) and i!=L*L-1:
        return   -Latticetemp[i]*(Latticetemp[i-1]+Latticetemp[i+1]+Latticetemp[i-L]+Latticetemp[i%L])
    if  i%L==(L-1) and i!=(L-1) and i!=L*L-1:
        return  -Latticetemp[i]*(Latticetemp[i-1]+Latticetemp[i+L]+Latticetemp[i-L]+Latticetemp[i-(L-1)])
    else:
        return  -Latticetemp[i]*(Latticetemp[i-1]+Latticetemp[i+1]+Latticetemp[i-L]+Latticetemp[i+L])

#Energy
Energy=0
for i in range(L*L):
    Energy=Energy+LocalHamiltonian(i,Lattice)

Energy=Energy/(2*L*L)
#print(Energy)

#Metropolis
def Metropolis(i,LatticeTemp1):
    LatticeTemp2=LatticeTemp1.copy()
    LatticeTemp2[i]=-LatticeTemp2[i]
    DeltaH=LocalHamiltonian(i,LatticeTemp2)-LocalHamiltonian(i,LatticeTemp1)
    if DeltaH<0:
        #print("accepted at step 1")
        return  LatticeTemp2
    else:
        if r.random()<mt.exp(-beta*(DeltaH)):
            #print("accepted at step 2 with prob=",mt.exp(-beta*(DeltaH)))
            return LatticeTemp2
        else:
            #print("non-accepted")
            return LatticeTemp1



#Magnetization
def Magnetization(Lattice):
    return np.mean(Lattice)


#Evolution
Betas=[]
MvsBeta=[]
for m in range(19):
    M=[]
    beta=(m+1)/(20)
    Betas.append(beta)
    #Termalization
    for j in range(1000*L*L):
        Lattice=Metropolis(r.randint(0,L**2-1),Lattice)
    for k in range(1000):
        for l in range(L*L):
            Lattice=Metropolis(r.randint(0,L**2-1),Lattice)
        M.append(abs(Magnetization(Lattice)))
    print("beta=",beta,"M=",np.mean(M))
    MvsBeta.append(abs(np.mean(M)))

plt.figure(1)
plt.plot(Betas,MvsBeta,"rd")
plt.show()
plt.close()

