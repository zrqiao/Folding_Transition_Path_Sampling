import argparse, math, random, gzip, pickle, types
import numpy as np
from functools import partial
import multiprocessing

from collections import defaultdict

def rate_symm(fi, fj):
    return math.exp(-(fj - fi) / 2)

def rate_metropolis(fi, fj):
    return math.exp(-min(0, fj - fi))

def rate_matrix_1d(F, rate_fcn, x_F_min, x_F_max):
    states = sorted(x for x in F if x >= x_F_min and x <= x_F_max)
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

def sample_transition_paths_1d(rates, id, save_time_only=False):
    N = rates.shape[0]
    i = 0
    tp = [(i, 0.)]
    ntps = 0
    while True:
        if i == 0 or i == N - 1:
            if tp[0][0] != i:# This path is a transition path
                if id % 1 == 0: print('step: %d \n' % (id))
                if not save_time_only:
                    #pickle.dump(tp, stream)
                    return tp
                else:
                    #pickle.dump((len(tp), sum(t[1] for t in tp)), stream)
                    return (len(tp), sum(t[1] for t in tp))
                ntps += 1

                if ntps >= npaths:
                    break
            tp = [(i, 0.)]
        if i == 0:
            j = i + 1
            cumulative_rate = rates[i,j]
        elif i == N - 1:
            j = i - 1
            cumulative_rate = rates[i,j]
        else:
            cumulative_rate = rates[i,i-1] + rates[i,i+1]
            r = cumulative_rate * np.random.random()
            if r < rates[i,i-1]:
                j = i - 1
            else:
                j = i + 1
        i = j
        dt = np.log(1 / random.random()) / cumulative_rate
        tp.append((i, dt))
# def transition_path_length_distribution(stream, nbins=10):
#     minL = min(len(tp) for tp in tps)
#     maxL = max(len(tp) for tp in tps) + 1
#     dL = (maxL - minL) / nbins
#     hist = np.zeros(nbins)
#     for tp in tps:
#         hist[int((len(tp) - minL) / dL)] += 1
#     return {minL + dL * (i + 0.5) : hist[i] / np.sum(hist) for i in range(nbins)}

def transition_path_time_distribution(stream, nbins=10):
    times = []
    while True:
        try:
            tp = pickle.load(stream)
            if isinstance(tp,list):
                times.append(sum(t[1] for t in tp))
                #print(times[-1])
            else:
                times.append(tp[1]) # Time information only
        except (IOError, EOFError):
            break
    print("Loaded %d transition paths" % len(times))
    def bootstrap_fcn_err(fcn, x, nsamples=100):
        return np.std([fcn(np.random.choice(x, len(x))) for i in range(nsamples)])
    meantime = np.mean(times)
    print("average time:", meantime, bootstrap_fcn_err(np.mean, times))
    times = [t / meantime for t in times]
    # stddevtime = np.std(times)
    # times = [t / stddevtime for t in times]
    mint = 0. # min(times)
    maxt = 32. # max(times) + 1 change!
    print("time stddev / mean:", np.std(times), bootstrap_fcn_err(np.std, times))
    dt = (maxt - mint) / nbins
    hist = np.zeros(nbins)
    for i in range(len(times)):
        hist[int((times[i] - mint) / dt)] += 1
    return {mint + dt * (i + 0.5) : (hist[i] / (dt * np.sum(hist)), np.sqrt(hist[i]) / (dt * np.sum(hist))) for i in range(nbins)}

def transition_path_time_cdf(streamin, streamout):
    times = []
    while True:
        try:
            tp = pickle.load(streamin)
            if isinstance(tp,list):
                times.append(sum(t[1] for t in tp))
            else:
                times.append(tp[1]) # Time information only
        except (IOError, EOFError):
            break
    print("Loaded %d transition paths" % len(times))
    times.sort()
    streamout.write("%g %g\n" % (0, 0))
    for i in range(len(times)):
        streamout.write("%g %g\n" % (times[i], (i + 1) / len(times)))

def jump_size_distribution(stream, dts, ddt=1, step_size=1):
    steps = {dt : defaultdict(int) for dt in dts}
    msd = {dt : 0. for dt in dts}
    ntps = 0
    while True:
        try:
            tp = pickle.load(stream)
            ntps += 1
        except (IOError, EOFError):
            break
        cumulative_time = sum(t[1] for t in tp)
        x_t = np.zeros(int(cumulative_time / ddt) + 2)
        i, t = 0, 0.
        for k in range(int(cumulative_time / ddt) + 2):
            while True:
                if i < len(tp) - 1 and \
                   math.fabs(t + tp[i + 1][1] - k * ddt) < math.fabs(t - k * ddt):
                    i += 1
                    t += tp[i][1]
                else:
                    x_t[k] = tp[i][0]
                    break
        for dt in dts:
            drift = (tp[-1][0] - tp[0][0]) / (cumulative_time / dt)
            for t in range(int(dt / ddt), len(x_t)):
                dx = (x_t[t] - x_t[t - int(dt / ddt)]) - drift #Should be - drift?
                msd[dt] += dx**2
                if dx > 0:
                    steps[dt][int(dx / step_size + 0.5)] += 1
                else:
                    steps[dt][int(dx / step_size - 0.5)] += 1
        del tp, x_t
    print("Loaded %d transition paths" % ntps)
    norm = {dt : sum(steps[dt].values()) for dt in dts}
    return {dt : msd[dt] / norm[dt] for dt in dts}, \
        {dt : {i * step_size : steps[dt][i] / norm[dt] \
               for i in range(min(steps[dt]), max(steps[dt]) + 1)} for dt in dts}

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('landscape', type=str, help="path to input landscape file")
    parser.add_argument('--stored-paths', type=str, default='tps.p.gz', \
                        help="path to stored paths [tps.p.gz]")
    parser.add_argument('--time-only', action='store_true', help="save length/time only [False]")
    clargs = parser.parse_args()

    with open(clargs.landscape, 'r') as f:
        F = {np.float(line.split()[0]) : np.float(line.split()[1]) for line in f \
             if len(line) > 1 and line[0] != '#'}
    x_F_mid = (min(F) + max(F)) / 2
    x_F_min = min((z for z in F.items() if z[0] < x_F_mid), key=lambda z: (z[1], z[0]))[0]
    x_F_max = min((z for z in F.items() if z[0] > x_F_mid), key=lambda z: (z[1], -z[0]))[0]
    x_F_barrier = max((z for z in F.items() if z[0] > x_F_min and z[0] < x_F_max), \
                      key=lambda z: z[1])[0]

    print(min(F), x_F_min, x_F_barrier, x_F_max, max(F))
    print("barrier:", (2 * F[x_F_barrier] - F[x_F_min] - F[x_F_max]) / 2)

    T, states = rate_matrix_1d(F, rate_symm, x_F_min, x_F_max)
    print("nstates =", len(states))
    pi = stationary(np.transpose(T))
    q = committors(T)
    m = pi * q * (1. - q)
    kreaction = reaction_rate(pi, T, [0], [-1])
    print("log(reaction rate) =", math.log(kreaction))
    print("Writing tpt_1d.dat")
    with open('tpt_1d.dat', 'w') as f:
        f.write("# x q^+(x) p(x|TP) p(TP|x)\n")
        for i in range(len(q)):
            f.write("%g %g %g %g\n" % (states[i], q[i], m[i] / sum(m), 2. * m[i] / pi[i]))

    npaths = 50000

    print("Sampling %d transition paths and writing to %s..." % (npaths, clargs.stored_paths))
    with gzip.open(clargs.stored_paths, 'wb') as f:
        print('cores='+str(multiprocessing.cpu_count()))
        task=partial(sample_transition_paths_1d, T, save_time_only=clargs.time_only)
        pool = multiprocessing.Pool()  # creates a pool of process, controls worksers
        # the pool.map only accepts one iterable, so use the partial function
        # so that we only need to deal with one variable.
        #task(npaths)
        A=list(pool.map(task, np.arange(npaths))) # make our results with a map call

        for obj in A:
            pickle.dump(obj, f)

        pool.close()  # we are not adding any more processes

    # print("Writing simulated_tp_length_distribution.dat")
    # with gzip.open(clargs.stored_paths, 'rb') as f_tps, \
    #      open('simulated_tp_length_distribution.dat', 'w') as f:
    #     for L,p in transition_path_length_distribution(f_tps).items():
    #         f.write("%g %g\n" % (L, p))
    print("Writing simulated_tp_time_distribution.dat")
    with gzip.open(clargs.stored_paths, 'rb') as f_tps, \
         open('simulated_tp_time_distribution.dat', 'w') as f:
        for t,p in transition_path_time_distribution(f_tps, nbins=400).items():
            f.write("%g %g %g\n" % (t, p[0], p[1]))
    print("Writing simulated_tp_time_cdf.dat")
    with gzip.open(clargs.stored_paths, 'rb') as f_tps, \
         open('simulated_tp_time_cdf.dat', 'w') as f:
         transition_path_time_cdf(f_tps, f)

    step_size = 1
    with gzip.open(clargs.stored_paths, 'rb') as f_tps:
        msd, steps = jump_size_distribution(f_tps, [2.**p for p in range(8)], step_size=step_size)
    print("Writing simulated_tp_msd.dat")
    with open('simulated_tp_msd.dat', 'w') as f:
        for dt in msd:
            f.write("%g %g\n" % (dt, msd[dt]))
    print("Writing simulated_tp_step_size_distribution.dat")
    with open('simulated_tp_step_size_distribution.dat', 'w') as f:
        for dt in steps:
            for i in range(min(steps[dt]), max(steps[dt]) + 1):
                f.write("%g %g %g\n" % (math.log(dt) / math.log(2), i * step_size, steps[dt][i]))
            f.write("\n")
