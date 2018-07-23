import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit

def processing(path,ax):
    with open(path+'/backward_jumping_probability.dat','r') as f:
        backward_prob=np.loadtxt(f).transpose()
        ax.plot(backward_prob[0], backward_prob[1], label=path, ls='solid')
        ax.legend(loc='best', fontsize='large')
    return backward_prob

if __name__ == '__main__':

    plt.style.use('bmh')
    fig = plt.figure(figsize=(8, 5))

    colors = [plt.cm.jet(lt) for lt in range(0,8)]
    fig.add_axes()

    #mpl.rcParams['axes.color_cycle'] = colors
    mpl.rcParams['axes.titlesize']=20
    mpl.rcParams['axes.titleweight'] = 15
    backward_prob=np.array([])


    fig = plt.figure(figsize=(8, 5))
    fig.add_axes()


    ax = fig.add_subplot(111)
    #ax.set_title('Transition Path Transit Time Distribution', fontsize=20)
    ax.set_xlabel(r'$discrete states$', fontsize=12.5)
    ax.set_ylabel(r'$p(i->i-1|TP)$', fontsize=12.5)
    ax.set_yscale('linear')
    ax.set_ylim(0.4, 0.6)
    #ax.set_xlim(0, 4)

    path = 'Harmonic'
    processing(path,ax)

    path = 'Periodic'
    processing(path,ax)

    path = '2Barriers'
    #processing(path,ax)

    path = '3Barriers'
    processing(path,ax)

    path = '4Barriers'
    #processing(path,ax)

    path = '5Barriers'
    #processing(path,ax)

    fig.savefig('Backward_jump_probability.eps')

    plt.show()