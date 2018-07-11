import numpy as np
import matplotlib.pyplot as plt

global A
A=5.37

X=np.arange(-75,76,1)

X_test=np.arange(-35,35,1)

def Toy_Cos(x):
    return A/2*np.cos(2* np.pi* x/ 151)

def Toy_Cos_test(x):
    return A/2*np.cos(2* np.pi* x/ 71)

def HighFreq1(x):
    return A/2*(- 4/5 *np.cos(4* np.pi* x/ 111)+1/5)

def HighFreq2(x):
    return A / 2 * ( 4 / 5 * np.cos(6 * np.pi * x / 109) + 1 / 5)

def HighFreq3(x):
    return A / 2 * (- 4 / 5 * np.cos(8 * np.pi * x / 107) + 1 / 5)

def HighFreq4(x):
    return A / 2 * ( 4 / 5 * np.cos(10 * np.pi * x / 106.5) + 1 / 5)

if __name__ == '__main__':

    Truc=23
    HarmToy=Toy_Cos(X)
    HarmToy_test = Toy_Cos_test(X_test)
    TwoBarrier=Toy_Cos(X)
    TwoBarrier[Truc:-Truc]=HighFreq1(X[Truc:-Truc])
    ThreeBarrier = Toy_Cos(X)
    ThreeBarrier[Truc:-Truc] = HighFreq2(X[Truc:-Truc])
    FourBarrier = Toy_Cos(X)
    FourBarrier[Truc:-Truc] = HighFreq3(X[Truc:-Truc])
    FiveBarrier = Toy_Cos(X)
    FiveBarrier[Truc:-Truc] = HighFreq4(X[Truc:-Truc])


    with open('harmonic_test.dat', 'w') as f:
        np.savetxt(f,np.array([X_test,HarmToy_test]).transpose())

    with open('harmonic.dat', 'w') as f:
        np.savetxt(f,np.array([X,HarmToy]).transpose())

    with open('1barrier.dat', 'w') as f:
        np.savetxt(f,np.array([X,TwoBarrier]).transpose())

    with open('3barrier.dat', 'w') as f:
        np.savetxt(f,np.array([X,ThreeBarrier]).transpose())

    with open('4barrier.dat', 'w') as f:
        np.savetxt(f,np.array([X,TwoBarrier]).transpose())


    fig = plt.figure(figsize=(10, 5))
    fig.add_axes()
    ax=fig.add_subplot(111)
    ax.set_xlabel('discrete state')
    ax.set_ylabel('F/kT')
    ax.plot(X,HarmToy,ls='solid',color='black',marker='o',label='harmonic')
    ax.plot(X, TwoBarrier, ls='solid', color='b',marker='^',label='2Barriers')
    ax.plot(X, ThreeBarrier, ls='solid', color='g', marker='s',label='3Barriers')
    ax.plot(X, FourBarrier, ls='solid', color='r', marker='p',label='4Barriers')
    ax.plot(X, FiveBarrier, ls='solid', color='y', marker='h',label='5Barriers')
    ax.legend(loc='best')
    plt.show()