import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from collections import defaultdict


def processing(path):
    with open(path+'simulated_tp_msd.dat','r') as f1:
        sqmsd=np.zeros(8)
        for line in f1.readlines():
            linedata=list(map(np.float,line.split()))
            sqmsd[int(np.log2(linedata[0]))]=np.sqrt(linedata[1])

    with open(path+'simulated_tp_step_size_distribution.dat','r') as f:
        disp = {lt: defaultdict(np.float) for lt in range(0, 8)}
        for line in f.readlines():
            linedata = list(map(np.float, line.split()))
            if len(linedata)!=0:
                disp[int(linedata[0])][linedata[1]/sqmsd[int(linedata[0])]]=linedata[2]*sqmsd[int(linedata[0])]
                #print(line)
    return disp

if __name__ == '__main__':

    #print (plt.style.available)
    plt.style.use('bmh')
    fig = plt.figure(figsize=(12, 18))

    colors = [plt.cm.jet(lt) for lt in range(0,8)]
    fig.add_axes()

    #mpl.rcParams['axes.color_cycle'] = colors
    mpl.rcParams['axes.titlesize']=12
    mpl.rcParams['axes.titleweight'] = 15
    #mpl.rcParams['mathtext.default'] = 'regular'
    TPTDistribution=np.array([])

    path = 'Harmonic/'
    disp=processing(path)
    ax = fig.add_subplot(321)
    ax.set_title('Harmonic')
    ax.set_xlabel(r'$\Delta x/ \sqrt{\langle \Delta x (\Delta t)^2 \rangle}$',fontsize=12.5)
    ax.set_ylabel(r'$p(\Delta x| \Delta t)$',fontsize=12.5)
    ax.set_yscale('log')
    ax.set_ylim(1e-6,1)
    ax.set_xlim(-10,10)
    for lt in range(1,8):
        ax.plot(disp[lt].keys(), disp[lt].values(),marker='o', ls='solid',label=r'$log_2(k_0 \Delta t)=%g$'% lt)
        print()
    ax.legend(loc='best',fontsize='x-large')

    path = '2Barriers/'
    disp=processing(path)
    ax = fig.add_subplot(322)
    ax.set_title('2Barriers')
    ax.set_xlabel(r'$\Delta x/ \sqrt{\langle \Delta x (\Delta t)^2 \rangle}$',fontsize=12.5)
    ax.set_ylabel(r'$p(\Delta x| \Delta t)$',fontsize=12.5)
    ax.set_yscale('log')
    ax.set_ylim(1e-6,1)
    ax.set_xlim(-10,10)
    for lt in range(1,8):
        ax.plot(disp[lt].keys(), disp[lt].values(),marker='o', ls='solid',label=r'$log_2(k_0 \Delta t)=%g$'% lt)
    ax.legend(loc='best',fontsize='x-large')

    path = '3Barriers/'
    disp=processing(path)
    ax = fig.add_subplot(323)
    ax.set_title('3Barriers')
    ax.set_xlabel(r'$\Delta x/ \sqrt{\langle \Delta x (\Delta t)^2 \rangle}$',fontsize=12.5)
    ax.set_ylabel(r'$p(\Delta x| \Delta t)$',fontsize=12.5)
    ax.set_yscale('log')
    ax.set_ylim(1e-6,1)
    ax.set_xlim(-10,10)
    for lt in range(1,8):
        ax.plot(disp[lt].keys(), disp[lt].values(),marker='o', ls='solid',label=r'$log_2(k_0 \Delta t)=%g$'% lt)
    ax.legend(loc='best',fontsize='x-large')

    path = '4Barriers/'
    disp=processing(path)
    ax = fig.add_subplot(324)
    ax.set_title('4Barriers')
    ax.set_xlabel(r'$\Delta x/ \sqrt{\langle \Delta x (\Delta t)^2 \rangle}$',fontsize=12.5)
    ax.set_ylabel(r'$p(\Delta x| \Delta t)$',fontsize=12.5)
    ax.set_yscale('log')
    ax.set_ylim(1e-6,1)
    ax.set_xlim(-10,10)
    for lt in range(1,8):
        ax.plot(disp[lt].keys(), disp[lt].values(),marker='o', ls='solid',label=r'$log_2(k_0 \Delta t)=%g$'% lt)
    ax.legend(loc='best',fontsize='x-large')

    path = '5Barriers/'
    disp=processing(path)
    ax = fig.add_subplot(325)
    ax.set_title('5Barriers')
    ax.set_xlabel(r'$\Delta x/ \sqrt{\langle \Delta x (\Delta t)^2 \rangle}$',fontsize=12.5)
    ax.set_ylabel(r'$p(\Delta x| \Delta t)$',fontsize=12.5)
    ax.set_yscale('log')
    ax.set_ylim(1e-6,1)
    ax.set_xlim(-10,10)
    for lt in range(1,8):
        ax.plot(disp[lt].keys(), disp[lt].values(),marker='o', ls='solid',label=r'$log_2(k_0 \Delta t)=%g$'% lt)
    ax.legend(loc='best',fontsize='x-large')

    path = 'Periodic/'
    disp=processing(path)
    ax = fig.add_subplot(326)
    ax.set_title('Periodic')
    ax.set_xlabel(r'$\Delta x/ \sqrt{\langle \Delta x (\Delta t)^2 \rangle}$',fontsize=12.5)
    ax.set_ylabel(r'$p(\Delta x| \Delta t)$',fontsize=12.5)
    ax.set_yscale('log')
    ax.set_ylim(1e-6,1)
    ax.set_xlim(-10,10)
    for lt in range(1,8):
        ax.plot(disp[lt].keys(), disp[lt].values(),marker='o', ls='solid',label=r'$log_2(k_0 \Delta t)=%g$'% lt)
    ax.legend(loc='best',fontsize='x-large')

    fig.savefig('Finite_time_Displacement.eps')
    #plt.colorbar()
    plt.show()