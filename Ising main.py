import numpy as np  
import scipy as sp
import matplotlib.pyplot as plt
import random as r
import math as mt

#---------------------------------SETTINGS-----------------------------------
L=15 #size of the lattice
beta=0.2 #inverse temperature value, when not specified otherwise


Lattice=[] #initialization of the first lattice (saved in a list)
for i in range(L**2):
    if  r.random()<0.5:
        Lattice.append(1)   
    else:
        Lattice.append(-1)

#--------------------------------FUNCTIONS----------------------------------

def NewLattice(L):
    Lattice=[] #initialization of the lattice, saved in a list
    for i in range(L**2):
        if  r.random()<0.5:
            Lattice.append(1)   
        else:
            Lattice.append(-1)
    return Lattice

#Local Hamiltonian (H of the single lattice site)
def LocalHamiltonian(i,inputLattice):
    if  i==0:
        return  -inputLattice[i]*(inputLattice[i+1]+inputLattice[i+L]+inputLattice[i+L-1]+inputLattice[L*(L-1)] )
    if  i==(L-1):
        return  -inputLattice[i]*(inputLattice[i-1]+inputLattice[0]+inputLattice[i+L]+inputLattice[L*L-1])
    if  i==L*(L-1):
        return  -inputLattice[i]*(inputLattice[i+1]+inputLattice[L*L-1]+inputLattice[0]+inputLattice[L*(L-2)])
    if  i==L*L-1:
        return  -inputLattice[i]*(inputLattice[L*L-2]+inputLattice[L*(L-1)]+inputLattice[L-1]+ inputLattice[L*(L-1)-1])
    if  i%L==0 and i!=0 and i!=L*(L-1):
        return  -inputLattice[i]*(inputLattice[i+1]+inputLattice[i+L-1]+inputLattice[i-L]+inputLattice[i+L] )
    if  mt.floor(i/L)==0 and i!=(L-1) and i!=0:
        return  -inputLattice[i]*(inputLattice[i-1]+inputLattice[i+1]+inputLattice[i+L]+inputLattice[L*(L-1)+i%L])
    if  mt.floor(i/L)==(L-1) and i!=L*(L-1) and i!=L*L-1:
        return   -inputLattice[i]*(inputLattice[i-1]+inputLattice[i+1]+inputLattice[i-L]+inputLattice[i%L])
    if  i%L==(L-1) and i!=(L-1) and i!=L*L-1:
        return  -inputLattice[i]*(inputLattice[i-1]+inputLattice[i+L]+inputLattice[i-L]+inputLattice[i-(L-1)])
    else:
        return  -inputLattice[i]*(inputLattice[i-1]+inputLattice[i+1]+inputLattice[i-L]+inputLattice[i+L])

#Energy
def Energy(inputLattice):
    temp=[]
    for i in range(L*L):
        temp=temp+inputLattice(i,Lattice)
    return  temp/2


#Metropolis step
def Metropolis(i,inputLattice):
    tempLattice=inputLattice.copy()
    tempLattice[i]=-tempLattice[i]
    DeltaH=LocalHamiltonian(i,tempLattice)-LocalHamiltonian(i,inputLattice)
    if DeltaH<0:
        #print("accepted at step 1")
        return  tempLattice
    else:
        if r.random()<mt.exp(-beta*(DeltaH)):
            #print("accepted at step 2 with prob=",mt.exp(-beta*(DeltaH)))
            return tempLattice
        else:
            #print("non-accepted")
            return inputLattice


#Magnetization
def Magnetization(Lattice):
    return np.mean(Lattice)

#----------------------------------------------------FIRST SIMULATION ----------------------------------------------
Betas=[] # Ranging beta
MvsBeta=[]# Mean Magnetization M values for Betas values
Nterm=1000 # Termalization cycles, each with L*L Metropolis iterations at random-selected sites
Ndata=1000 # Number of data points collected, each single data point is collected after L*L Metropolis iterations at random-selected sites.
# Range of betas
BetaMin=1/10
BetaMax=10/10
Nsteps=20
M=[]#Magnetization
for m in range(Nsteps):
    Lattice=NewLattice(L)
    beta=BetaMin+m*((BetaMax-BetaMin)/Nsteps)
    Mtemp=[] # Magnetization at this beta 
    Betas.append(beta)
    for j in range(Nterm*L*L):  #Termalization cycles
        Lattice=Metropolis(r.randint(0,L**2-1),Lattice)
    for k in range(Ndata):      #Simulation cycles
        for l in range(L*L):
            Lattice=Metropolis(r.randint(0,L**2-1),Lattice)
        Mtemp.append(abs(Magnetization(Lattice)))
    print("beta=",beta,"M=",np.mean(Mtemp))
    M.append(Mtemp) #Saving all the data (lattice configurations) in a list(data point) of list(whole lattice)
    MvsBeta.append(abs(np.mean(Mtemp))) #(Just for a quick visualization of results)
 

plt.figure(1) #Very quick visualization of the results
plt.plot(Betas,MvsBeta,"r+")
plt.xlabel("Beta")
plt.ylabel("M")
plt.show()
plt.close()
