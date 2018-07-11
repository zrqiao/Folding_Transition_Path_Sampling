import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    TPTDistribution=np.array([])

    with open('simulated_tp_time_distribution.dat','r') as f:
        TPTDistribution=np.loadtxt(f).transpose()

    fig = plt.figure(figsize=(10, 5))
    fig.add_axes()
    ax = fig.add_subplot(111)
    ax.set_xlabel('time')
    ax.set_ylabel('P')
    ax.set_yscale('log')
    ax.plot(TPTDistribution[0], TPTDistribution[1], ls='solid', color='black', marker='o')
    ax.legend(loc='best')
    plt.show()