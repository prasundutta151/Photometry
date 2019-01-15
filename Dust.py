import numpy as np
import  matplotlib.pyplot as plt
import scipy as sp

CC = 3.e8

def Filter(lam, par):

     return 1.*np.exp(-(lam-par[0])**2/2./par[1]**2)


def BBT(lam, TT):

    return 1./(lam**5*(np.exp(1.43e7/TT/lam)-1))


def DustA(lam):

     BBV =  lam**4*np.exp(-lam/30)
     BBV = BBV/BBV.max()
     BBV = BBV + 10000*np.exp(-((lam-80)/8)**2)
     BBV = BBV + 30*np.exp(-((lam-220)/20)**2)
     BBV = np.log(BBV)
     BBV = 18*BBV/np.max(BBV) + 50
     BBV = BBV/70.*12
     return BBV


lam = np.linspace(10, 1000., 1024)
lamU0 = 365.
sigU  = 68.
lamB0 = 440.
sigB  = 98.
lamV0 = 550.
sigV  = 89.

parU  = [lamU0, sigU]
parB  = [lamB0, sigB]
parV  = [lamV0, sigV]

plt.ion()


RU = Filter(lam, parU)
RB = Filter(lam, parB)
RV = Filter(lam, parV)

Temp = np.linspace(3000, 20000, 128)
CBV  = np.zeros(128)
CUB  = np.zeros(128)
CBVD  = np.zeros(128)
CUBD  = np.zeros(128)

Av = np.exp(-DustA(lam))     

for ii in range(128):

     BB   = BBT(lam, Temp[ii])
     BB = BB/np.max(BB)
     
  #   plt.plot(lam, RU, color='blue', linewidth=2)
  #   plt.plot(lam, RB, color='green', linewidth=2)
  #   plt.plot(lam, RV, color='red', linewidth=2)
  #   plt.plot(lam, BB, color='black', linewidth=4)
     
     
     MU = np.sum(RU*BB)
     MB = np.sum(RB*BB)
     MV = np.sum(RV*BB)


     MUD = np.sum(RU*BB*Av)
     MBD = np.sum(RB*BB*Av)
     MVD = np.sum(RV*BB*Av)

     CBmV = -2.5*np.log10(MU/MV)
     CUmB = -2.5*np.log10(MU/MB)
   
     CBmVD = -2.5*np.log10(MUD/MVD)
     CUmBD = -2.5*np.log10(MUD/MBD)
     
     CBV[ii] = CBmV
     CUB[ii] = CUmB

     CBVD[ii] = CBmVD
     CUBD[ii] = CUmBD



plt.clf()
plt.xlim(-1.4, 1.5)
plt.ylim(1.2, -0.2)
plt.plot(CBV,CUB)
plt.plot(CBVD, CUBD)
plt.show()

