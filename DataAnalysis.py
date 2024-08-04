import numpy as np  
import scipy as sp
import matplotlib.pyplot as plt
import random as r
import math as mt

M=np.loadtxt("data1.gz", dtype='float')
Betas=np.loadtxt("Betas1.gz", dtype='float')


def Jacknife(Vec):
    lenght=np.size(Vec)
    VecCopies=[]
    for j in range(lenght):
        temp=np.delete(Vec,r.randint(0,lenght-1))
        VecCopies.append(temp)
    #Means=np.mean(VecCopies,axis=0)
    Means = np.array([np.mean(x) for x in VecCopies])
    return [np.mean(Means),np.std(Means)]

def Bootstrap(Vec):
    lenght=np.size(Vec)
    VecCopies=[]
    for j in range(lenght):
        indices=np.random.randint(0,lenght,lenght)
        VecCopies.append(Vec[0,indices])
    Means = np.array([np.mean(x) for x in VecCopies])
    return [np.mean(Means),np.std(Means)]


def Bootstrap2(Vec,Delta):
    lenght=np.size(Vec)
    N=lenght//Delta
    VecCopies=[]
    print(Delta)
    for i in range(lenght):
        indices=np.random.randint(0,lenght-Delta,N)
        indice = np.array([ x+j for x in indices for j in range(Delta)])

        VecCopies.append(Vec[0,indice])
    Means = np.array([np.mean(x) for x in VecCopies])
    return [np.mean(Means),np.std(Means)]                
        


results=[]
for i in range(np.size(M,axis=0)):
    results.append(Bootstrap2(np.take(M,[i],axis=0),100))
print(results)
MeanM=np.take(results,[0],axis=1)
MeanM=np.reshape(MeanM,np.size(MeanM))
DeltaM=np.take(results,[1],axis=1)
DeltaM=np.reshape(DeltaM,np.size(DeltaM))
print(MeanM)
print(DeltaM)

plt.figure(1) 
plt.errorbar(Betas,MeanM,yerr=DeltaM,marker =".")
plt.xlabel("Beta")
plt.ylabel("MeanM")
plt.show()
plt.close()
