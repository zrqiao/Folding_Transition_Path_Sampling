import numpy as np

global A
A=5.37

X=np.arange(-35,35,1)

def Toy_Cos(x):
    return A*np.cos(2* np.pi* x/ 71)

if __name__ == '__main__':

    HarmToy=Toy_Cos(X)
    with open('harmonic.dat', 'w') as f:
        np.savetxt(f,np.array([X,HarmToy]).transpose())