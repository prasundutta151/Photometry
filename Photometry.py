import numpy as np
import  matplotlib.pyplot as plt
import scipy as sp

CC = 3.e8

def Filter(lam, par):

     return 1.*np.exp(-(lam-par[0])**2/2./par[1]**2)


def BBT(lam, TT):

    return 1./(lam**5*(np.exp(1.43e7/TT/lam)-1))



lam = np.linspace(10, 1600., 1024)
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
plt.figure()

RU = Filter(lam, parU)
RB = Filter(lam, parB)
RV = Filter(lam, parV)

Temp = np.linspace(4000, 20000, 128)
CBV  = np.zeros(128)
CUB  = np.zeros(128)

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
     
     CBmV = -2.5*np.log10(MU/MV)
     CUmB = -2.5*np.log10(MU/MB)
     
     CBV[ii] = CBmV
     CUB[ii] = CUmB


plt.clf()
plt.xlim(-1.5, 1.5)
plt.ylim(1.2, -0.2)
plt.plot(CBV, CUB)
plt.show()

