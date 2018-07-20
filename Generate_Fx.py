import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

global A
A=5.37
Truc = 23
X=np.arange(-75,76,1)
X_l=np.arange(-75,75+0.1,0.1)

X_test=np.arange(-35,35,1)

def Toy_Cos(x):
    return A/2*np.cos(2* np.pi* x/ 151)

def Toy_Cos_flat(x):
    return A/8*np.cos(2* np.pi* x/ 151)

def Toy_Periodic(x):
    return (A/2*np.cos(2* np.pi* x/ 151)+2*np.cos(10*np.pi* x/151))*(A/(A+4))

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

def HighFreq5(x):
    return A / 2 * ( 4 / 5 * np.cos(100 * np.pi * x / 106.5) + 1 / 5)



if __name__ == '__main__':


    HarmToy=Toy_Cos(X)
    Harm_l=Toy_Cos(X_l)
    Harm_flat = Toy_Cos_flat(X)
    Harm_flat_l = Toy_Cos_flat(X_l)
    HarmToy_test = Toy_Cos_test(X_test)
    TwoBarrier=Toy_Cos(X)
    TwoBarrier[Truc:-Truc]=HighFreq1(X[Truc:-Truc])
    TwoBarrier_l = Toy_Cos(X_l)
    TwoBarrier_l[Truc*10:-Truc*10] = HighFreq1(X_l[Truc*10:-Truc*10])
    ThreeBarrier = Toy_Cos(X)
    ThreeBarrier[Truc:-Truc] = HighFreq2(X[Truc:-Truc])
    ThreeBarrier_l = Toy_Cos(X_l)
    ThreeBarrier_l[Truc*10:-Truc*10] = HighFreq2(X_l[Truc*10:-Truc*10])
    FourBarrier = Toy_Cos(X)
    FourBarrier[Truc:-Truc] = HighFreq3(X[Truc:-Truc])
    FourBarrier_l = Toy_Cos(X_l)
    FourBarrier_l[Truc*10:-Truc*10] = HighFreq3(X_l[Truc*10:-Truc*10])
    FiveBarrier = Toy_Cos(X)
    FiveBarrier[Truc:-Truc] = HighFreq4(X[Truc:-Truc])
    FiveBarrier_l = Toy_Cos(X_l)
    FiveBarrier_l[Truc*10:-Truc*10] = HighFreq4(X_l[Truc*10:-Truc*10])
    Highfreq_l = Toy_Cos(X_l)
    Highfreq_l[Truc*10:-Truc*10] = HighFreq5(X_l[Truc*10:-Truc*10])
    Periodic=Toy_Periodic(X)
    Periodic_l = Toy_Periodic(X_l)

    with open('harmonic_test/PES.dat', 'w') as f:
        np.savetxt(f,np.array([X_test,HarmToy_test]).transpose())

    with open('Harmonic/PES_l.dat', 'w') as f:
        np.savetxt(f,np.array([X_l,Harm_l]).transpose())

    with open('Harmonic_flat/PES.dat', 'w') as f:
        np.savetxt(f,np.array([X,Harm_flat]).transpose())

    with open('Harmonic_flat/PES_l.dat', 'w') as f:
        np.savetxt(f,np.array([X_l,Harm_flat_l]).transpose())

    with open('2Barriers/PES_l.dat', 'w') as f:
        np.savetxt(f,np.array([X_l,TwoBarrier_l]).transpose())

    with open('3Barriers/PES_l.dat', 'w') as f:
        np.savetxt(f,np.array([X_l,ThreeBarrier_l]).transpose())

    with open('4Barriers/PES_l.dat', 'w') as f:
        np.savetxt(f,np.array([X_l,FourBarrier_l]).transpose())

    with open('5Barriers/PES_l.dat', 'w') as f:
        np.savetxt(f,np.array([X_l,FiveBarrier_l]).transpose())

    with open('Periodic/PES_l.dat', 'w') as f:
        np.savetxt(f,np.array([X_l,Periodic_l]).transpose())

    with open('Highfreq/PES_l.dat', 'w') as f:
        np.savetxt(f,np.array([X_l,Highfreq_l]).transpose())

    with open('Harmonic/PES.dat', 'w') as f:
        np.savetxt(f,np.array([X,HarmToy]).transpose())

    with open('2Barriers/PES.dat', 'w') as f:
        np.savetxt(f,np.array([X,TwoBarrier]).transpose())

    with open('3Barriers/PES.dat', 'w') as f:
        np.savetxt(f,np.array([X,ThreeBarrier]).transpose())

    with open('4Barriers/PES.dat', 'w') as f:
        np.savetxt(f,np.array([X,FourBarrier]).transpose())

    with open('5Barriers/PES.dat', 'w') as f:
        np.savetxt(f,np.array([X,FiveBarrier]).transpose())

    with open('Periodic/PES.dat', 'w') as f:
        np.savetxt(f,np.array([X,Periodic]).transpose())

    mpl.rcParams['axes.titlesize']=40
    mpl.rcParams['axes.titleweight'] = 15
    plt.style.use('ggplot')
    fig = plt.figure(figsize=(16, 10))

    fig.add_axes()
    ax=fig.add_subplot(231)
    ax.set_title('Free Energy Landscape')
    ax.set_xlabel('discrete state',fontsize='15')
    ax.set_ylabel(r'$F/kT$',fontsize='15')
    ax.plot(X,HarmToy,ls='solid',marker='o',label='Harmonic')
    ax.plot(X, TwoBarrier, ls='solid',marker='^',label='2Barriers')
    ax.legend(loc='best',fontsize='xx-large')
    ax=fig.add_subplot(232)
    ax.set_title('Free Energy Landscape')
    ax.set_xlabel('discrete state',fontsize='15')
    ax.set_ylabel(r'$F/kT$',fontsize='15')
    ax.plot(X,HarmToy,ls='solid',marker='o',label='Harmonic')
    ax.plot(X, ThreeBarrier, ls='solid', marker='s',label='3Barriers')
    ax.legend(loc='best',fontsize='xx-large')
    ax=fig.add_subplot(233)
    ax.set_title('Free Energy Landscape')
    ax.set_xlabel('discrete state',fontsize='15')
    ax.set_ylabel(r'$F/kT$',fontsize='15')
    ax.plot(X,HarmToy,ls='solid',marker='o',label='Harmonic')
    ax.plot(X, FourBarrier, ls='solid', marker='p',label='4Barriers')
    ax.legend(loc='best',fontsize='xx-large')
    ax=fig.add_subplot(234)
    ax.set_title('Free Energy Landscape')
    ax.set_xlabel('discrete state',fontsize='15')
    ax.set_ylabel(r'$F/kT$',fontsize='15')
    ax.plot(X,HarmToy,ls='solid',marker='o',label='Harmonic')
    ax.plot(X, FiveBarrier, ls='solid', marker='h',label='5Barriers')
    ax.legend(loc='best',fontsize='xx-large')
    ax=fig.add_subplot(235)
    ax.set_title('Free Energy Landscape')
    ax.set_xlabel('discrete state',fontsize='15')
    ax.set_ylabel(r'$F/kT$',fontsize='15')
    ax.plot(X,HarmToy,ls='solid',marker='o',label='Harmonic')
    ax.plot(X, Periodic, ls='solid', marker='h',label='Periodic')
    ax.legend(loc='best',fontsize='xx-large')

    fig. savefig('FELs.eps')
    plt.show()

    fig2 = plt.figure(figsize=(8, 5))

    ax=fig2.add_subplot(111)
    ax.set_title('Free Energy Landscape')
    ax.set_xlabel('discrete state',fontsize='15')
    ax.set_ylabel(r'$F/kT$',fontsize='15')
    ax.plot(X_l,Harm_l,ls='solid',label='Harmonic')
    ax.plot(X_l,Harm_flat_l,ls='solid',label='Harmonic_flat')
    ax.plot(X_l,Highfreq_l, ls='solid',label='HighFreq')
    ax.legend(loc='best',fontsize='xx-large')
    fig2. savefig('FELs_2.eps')


    plt.show()