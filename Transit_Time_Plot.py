import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit

def Theo_Distribution(t,F,omega):
    return omega*F*np.exp(-(F/(np.sinh(omega*t/2))**2))*np.cosh(omega*t/2)/(np.sinh(omega*t/2))**3

def Theo_Tail(t,omega,m):
    return m*np.exp(-omega*t)

def processing(path,ax):
    with open(path+'/simulated_tp_time_distribution.dat','r') as f:
        TPTDistribution=np.loadtxt(f).transpose()
        paras0,cov0=curve_fit(Theo_Tail,TPTDistribution[0][20:],TPTDistribution[1][20:])
        omega=paras0[0]
        domega=np.sqrt(np.diag(cov0))[0]
        paras,cov=curve_fit(lambda t, F:Theo_Distribution(t,F,omega),TPTDistribution[0],TPTDistribution[1],maxfev=50000)
        F=paras[0]
        dF=np.sqrt(np.diag(cov))[0]
        ax.plot(TPTDistribution[0], TPTDistribution[1], label=path, ls='solid', marker='o')
        ax.plot(TPTDistribution[0],Theo_Distribution(TPTDistribution[0],F,omega),ls='--',c='black',label=r'$\omega =%5.2f $, $\Delta F=%5.2f kT $' % (omega,F))
        ax.legend(loc='best', fontsize='large')
    return TPTDistribution

if __name__ == '__main__':

    plt.style.use('ggplot')
    fig = plt.figure(figsize=(8, 5))

    colors = [plt.cm.jet(lt) for lt in range(0,8)]
    fig.add_axes()

    #mpl.rcParams['axes.color_cycle'] = colors
    mpl.rcParams['axes.titlesize']=20
    mpl.rcParams['axes.titleweight'] = 15
    TPTDistribution=np.array([])


    fig = plt.figure(figsize=(10, 7.5))
    fig.add_axes()


    ax = fig.add_subplot(111)
    #ax.set_title('Transition Path Transit Time Distribution', fontsize=20)
    ax.set_xlabel(r'$t_{AB}/\langle t_{AB} \rangle$', fontsize=12.5)
    ax.set_ylabel(r'$p(t_{AB})$', fontsize=12.5)
    ax.set_yscale('log')
    ax.set_ylim(5e-3, 2)
    ax.set_xlim(0, 4)

    path = 'Harmonic'
    processing(path,ax)

    path = 'Periodic'
    processing(path,ax)

    path = '2Barriers'
    processing(path,ax)

    path = '3Barriers'
    processing(path,ax)

    path = '4Barriers'
    processing(path,ax)

    path = '5Barriers'
    processing(path,ax)


    fig.savefig('Transit_Time_distribution_log.eps')


    plt.show()