import argparse, math, random, gzip, pickle, types
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from collections import defaultdict

def rate_symm(fi, fj):
    return math.exp(-(fj - fi) / 2)

def gibbs(fi, fj):
    return np.exp(-(fj - fi))

def rate_matrix_1d(F, rate_fcn, x_F_min, x_F_max):
    states = sorted(x for x in F)
    N = len(states)
    rates = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if abs(j - i) == 1:
                rates[i,j] = rate_fcn(F[states[i]], F[states[j]])
    for i in range(N):
        rates[i,i] = -sum(rates[i,:])
    return rates, states

def committors(rates, source=[0], sink=[-1], reverse=False):
    if 0 not in source or len(source) != max(source) + 1 or -1 not in sink or len(sink) != -min(sink):#Format
        raise Exception("invalid order of source and sink states")
    N = rates.shape[0]
    q = np.zeros(N)#States number
    if not reverse:
        b = np.array([-sum(rates[i,j] for j in sink) for i in range(len(source),N-len(sink))])#Check
        for j in sink:
            q[j] = 1.
    else:
        # Check the following line:
        b = np.array([-sum(rates[i,j] for j in source) for i in range(len(source),N-len(sink))])
        for j in source:
            q[j] = 1.
    q[max(source)+1:min(sink)] \
        = np.linalg.solve(rates[max(source)+1:min(sink),max(source)+1:min(sink)], b)
    return q

def stationary(rates):
    eig = np.linalg.eig(rates)
    pi = np.fabs(eig[1][:,np.argmax(eig[0])])
    return pi / np.sum(pi)

def reaction_rate(pi, T, source, sink):
    q = committors(T, source=source, sink=sink)
    k = sum(sum(pi[i] * T[i,j] * (q[j] - q[i]) \
                for j in range(len(pi)) if j not in source and q[j] > q[i]) \
            for i in range(len(pi)) if i in source)
    return k

def processing_reactive_prob(path, ax):

    with open(path+'/PES_l.dat', 'r') as f:
        F = {float(line.split()[0]) : float(line.split()[1]) for line in f \
             if len(line) > 1 and line[0] != '#'}
    x_F_mid = (min(F) + max(F)) / 2
    x_F_min = min((z for z in F.items() if z[0] < x_F_mid), key=lambda z: (z[1], z[0]))[0]
    x_F_max = min((z for z in F.items() if z[0] > x_F_mid), key=lambda z: (z[1], -z[0]))[0]
    x_F_barrier = max((z for z in F.items() if z[0] > x_F_min and z[0] < x_F_max), \
                      key=lambda z: z[1])[0]

    state_xs=np.array(list(F.keys()))
    state_Fs=np.array(list(F.values()))
    ref=np.zeros(len(state_Fs))
    #PDF=gibbs(ref,state_Fs)#Probability Distribution Function
    #Z=np.sum(PDF)#Partition Function
    #PDF=PDF/Z

    T, states = rate_matrix_1d(F, rate_symm, x_F_min, x_F_max)
    print("nstates =", len(states))
    pi = stationary(np.transpose(T))
    q = committors(T, [0], [-1])
    m = pi * q * (1. - q)
    #kreaction = reaction_rate(pi, T, [i for i in range (460)], [-i for i in range (1,461)])
    #print("log(reaction rate) =", math.log(kreaction))

    PDF=m/np.sum(m)*len(m) #q distribution by x
    tp_x=2. * m / pi
    dq=0.025

    p_q = np.zeros(int(1/dq)+1)
    q_coor=np.arange(0.,1.+dq,dq)

    for i in range(1,len(q)-1):
        p_q[int(np.round(q[i] / dq))] +=PDF[i]*dq

    with open('TP_distribution_'+path+'dat', 'w') as f:
       np.savetxt(f,p_q)
    ax.plot(state_xs, tp_x, label=path, ls='solid')

def processing_q(path, ax):

    with open(path+'/PES.dat', 'r') as f:
        F = {float(line.split()[0]) : float(line.split()[1]) for line in f \
             if len(line) > 1 and line[0] != '#'}
    x_F_mid = (min(F) + max(F)) / 2
    x_F_min = min((z for z in F.items() if z[0] < x_F_mid), key=lambda z: (z[1], z[0]))[0]
    x_F_max = min((z for z in F.items() if z[0] > x_F_mid), key=lambda z: (z[1], -z[0]))[0]
    x_F_barrier = max((z for z in F.items() if z[0] > x_F_min and z[0] < x_F_max), \
                      key=lambda z: z[1])[0]

    state_xs=np.array(list(F.keys()))
    state_Fs=np.array(list(F.values()))
    ref=np.zeros(len(state_Fs))
    #PDF=gibbs(ref,state_Fs)#Probability Distribution Function
    #Z=np.sum(PDF)#Partition Function
    #PDF=PDF/Z

    T, states = rate_matrix_1d(F, rate_symm, x_F_min, x_F_max)
    print("nstates =", len(states))
    pi = stationary(np.transpose(T))
    q = committors(T, [0], [-1])
    m = pi * q * (1. - q)
    #kreaction = reaction_rate(pi, T, [i for i in range (460)], [-i for i in range (1,461)])
    #print("log(reaction rate) =", math.log(kreaction))

    PDF=m/np.sum(m)*len(m) #q distribution by x
    tp_x=2. * m / pi
    dq=0.025

    p_q = np.zeros(int(1/dq)+1)
    q_coor=np.arange(0.,1.+dq,dq)
    recrossing=np.zeros(len(state_xs))


    for i in range(1,len(q)-1):
        p_q[int(np.round(q[i] / dq))] +=PDF[i]*dq
    for i in range(len(state_xs)-1):
        recrossing[i]=T[i,i-1]*q[i-1]/(T[i,i-1]*q[i-1]+T[i,i+1]*q[i+1])

    with open('recrossing_distribution_'+path+'dat', 'w') as f:
       np.savetxt(f,q)
    ax.plot(state_xs, recrossing, label=path, ls='solid')

def processing_recrossing(path,ax):

    with open(path+'/PES_l.dat', 'r') as f:
        F = {float(line.split()[0]) : float(line.split()[1]) for line in f \
             if len(line) > 1 and line[0] != '#'}
    x_F_mid = (min(F) + max(F)) / 2
    x_F_min = min((z for z in F.items() if z[0] < x_F_mid), key=lambda z: (z[1], z[0]))[0]
    x_F_max = min((z for z in F.items() if z[0] > x_F_mid), key=lambda z: (z[1], -z[0]))[0]
    x_F_barrier = max((z for z in F.items() if z[0] > x_F_min and z[0] < x_F_max), \
                      key=lambda z: z[1])[0]

    state_xs=np.array(list(F.keys()))
    state_Fs=np.array(list(F.values()))
    ref=np.zeros(len(state_Fs))
    #PDF=gibbs(ref,state_Fs)#Probability Distribution Function
    #Z=np.sum(PDF)#Partition Function
    #PDF=PDF/Z

    T, states = rate_matrix_1d(F, rate_symm, x_F_min, x_F_max)
    print("nstates =", len(states))
    pi = stationary(np.transpose(T))
    q = committors(T, [0], [-1])
    m = pi * q * (1. - q)
    #kreaction = reaction_rate(pi, T, [i for i in range (460)], [-i for i in range (1,461)])
    #print("log(reaction rate) =", math.log(kreaction))

    PDF=m/np.sum(m)*len(m) #q distribution by x
    tp_x=2. * m / pi
    dq=0.025

    p_q = np.zeros(int(1/dq)+1)
    q_coor=np.arange(0.,1.+dq,dq)

    for i in range(1,len(q)-1):
        p_q[int(np.round(q[i] / dq))] +=PDF[i]*dq

    with open('q_distribution_'+path+'dat', 'w') as f:
       np.savetxt(f,q)
    ax.plot(state_xs, q, label=path, ls='solid')

if __name__ == '__main__':


    plt.style.use('bmh')

    fig = plt.figure(figsize=(16, 12))

    fig.add_axes()

    #mpl.rcParams['axes.color_cycle'] = colors
    mpl.rcParams['axes.titlesize']=20
    mpl.rcParams['axes.titleweight'] = 15
    TPTDistribution=np.array([])


    fig = plt.figure(figsize=(8, 5))
    fig.add_axes()

    #ax = fig.add_subplot(111)
    #ax.set_xlabel('discrete states', fontsize=12.5)
    #ax.set_ylabel(r'$p(TP|x_{state})$', fontsize=12.5)
    #ax.set_yscale('linear')
    #ax.set_ylim(5e-3, 2)



    #processing('Harmonic',ax)
    #processing('Harmonic_flat',ax)
    #processing('Periodic',ax)
    #processing('2Barriers',ax)
    #processing('3Barriers',ax)
    #processing('4Barriers',ax)
    #processing('5Barriers',ax)
    #processing('Highfreq',ax)

    #ax.plot(state_xs,PDF1,state_xs,PDF2)
    #ax.plot(state_xs,q1,state_xs,q2)
    #ax.legend(loc='best', fontsize='x-large')
    #plt.show()
    #fig.savefig('TP_Distribution.png')

    fig1 = plt.figure(figsize=(8, 5))
    fig1.add_axes()

    ax1 = fig1.add_subplot(111)
    ax1.set_xlabel('discrete states', fontsize=12.5)
    ax1.set_ylabel(r'$P(i->i-1|TP)$', fontsize=12.5)
    ax1.set_yscale('linear')
    ax1.set_ylim(0.25, .75)

    processing_q('Harmonic', ax1)
    #processing_q('Harmonic_flat', ax1)
    processing_q('Periodic',ax1)
    processing_q('2Barriers',ax1)
    processing_q('3Barriers',ax1)
    processing_q('4Barriers',ax1)
    processing_q('5Barriers',ax1)
    #processing_q('Highfreq', ax1)

    #ax.plot(state_xs,PDF1,state_xs,PDF2)
    #ax.plot(state_xs,q1,state_xs,q2)
    ax1.legend(loc='best', fontsize='large')
    plt.show()
    fig.savefig('Recrossing_Distribution.png')


'''
    path = '2Barriers/'
    TPTDistribution=processing(path)
    ax.plot(TPTDistribution[0], TPTDistribution[1], label='2Barriers', ls='solid', marker='o')
    ax.legend(loc='best', fontsize='x-large')

    path = '3Barriers/'
    TPTDistribution=processing(path)
    ax.plot(TPTDistribution[0], TPTDistribution[1], label='3Barriers', ls='solid', marker='o')
    ax.legend(loc='best', fontsize='x-large')

    path = '4Barriers/'
    TPTDistribution=processing(path)
    ax.plot(TPTDistribution[0], TPTDistribution[1], label='4Barriers', ls='solid', marker='o')
    ax.legend(loc='best', fontsize='x-large')

    path = '5Barriers/'
    TPTDistribution=processing(path)
    ax.plot(TPTDistribution[0], TPTDistribution[1], label='5Barriers', ls='solid', marker='o')
    ax.legend(loc='best', fontsize='x-large')

    fig.savefig('Transit_Time_distirbution_linear.png')

'''
