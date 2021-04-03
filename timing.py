from orbit3d import orbit
import numpy as np
import time


def calc_times(n, ecc, N=None, n_rand=100):

    tstart = time.time()
    if N is None:
        N = 1000000//n
    
    t0 = 2455197.5
    
    alldata = []
    for k in range(n_rand):
        data = orbit.Data(1, 'HGCA_vDR2_corrected.fits', '', '', verbose=False)
        data.custom_epochs(np.sort(np.random.rand(n))*365.25 + t0, iplanet=0)
        data.nTot = n
        alldata += [data]
    
    model = orbit.Model(alldata[0])

    params_list = [1, 1e-5, 1, 0, ecc**0.5, 0, 1, 1, 1]
    par = orbit.Params(params_list)
    
    t1 = time.perf_counter_ns()
    for i in range(N//n_rand):
        for k in range(n_rand):
            orbit.calc_EA_RPP(alldata[k], par, model)
    t_RPP = (time.perf_counter_ns() - t1)/(N*n)
    
    t1 = time.perf_counter_ns()
    Nit = 0
    for i in range(N//n_rand):
        for k in range(n_rand):
            Nit = orbit.calc_EA_compare(alldata[k], par, model, Nit=Nit, solver='g')
    t_goatherd = (time.perf_counter_ns() - t1)/(N*n)
    t1 = time.perf_counter_ns()
    for i in range(N//n_rand):
        for k in range(n_rand):
            orbit.calc_EA_compare(alldata[k], par, model, solver='r')
    t_radvel = (time.perf_counter_ns() - t1)/(N*n)

    t1 = time.perf_counter_ns()
    for i in range(N//n_rand):
        for k in range(n_rand):
            orbit.calc_EA_compare(alldata[k], par, model, solver='b')
            
    t_batman = (time.perf_counter_ns() - t1)/(N*n)

    t1 = time.perf_counter_ns()
    for i in range(N//n_rand):
        for k in range(n_rand):
            orbit.calc_EA_compare(alldata[k], par, model, solver='e')
            
    t_nijenhuis = (time.perf_counter_ns() - t1)/(N*n)

    t1 = time.perf_counter_ns()
    for k in range(N):
        orbit.sincos(model, alldata[0].nTot)
    t_sincos = (time.perf_counter_ns() - t1)/(N*n)

    return [t_RPP, t_goatherd, t_radvel, t_batman, t_nijenhuis, t_sincos]


if __name__ == "__main__":

    ecc = np.asarray([1e-4, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
                      0.75, 0.8, 0.85, 0.9, 0.95, 0.98, 0.99])

    n = np.asarray([25, 100, 500])
    
    for _n in n:
        for _ecc in ecc:
            tvals = []
            for k in range(3):
                tvals += [calc_times(_n, _ecc, N=1000000//_n, n_rand=50)]
            print('%3d %6.2g %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f' %
                  tuple([_n, _ecc] + list(np.amin(np.asarray(tvals), axis=0))))
