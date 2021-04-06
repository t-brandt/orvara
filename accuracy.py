from orbit3d import orbit
import numpy as np
import time

def calc_accuracy(n, ecc):

    data = orbit.Data(0, '', '', '', verbose=False)

    # Define the right answers and convert to mean anomaly
    
    EA_in = (2*np.random.rand(n) - 1)*np.pi
    sEA_in = np.sin(EA_in)
    cEA_in = np.cos(EA_in)
    M = EA_in - ecc*sEA_in

    data.custom_epochs(M, refep=0, iplanet=0)
    data.nTot = n    
    model = orbit.Model(data)

    # Hijack the parameters class: desired eccentricity, period=2pi.
    par = orbit.Params([1]*9)
    par.ecc = ecc
    par.per = 2*np.pi
    par.lam = par.arg = 0

    # The last argument prevents roundoff error from incrementing
    # mean anomaly.  For a large number of epochs this eps*sqrt(n)
    # accumulation dominates error from the Kepler solver.
    
    orbit.calc_EA_RPP(data, par, model, True)
    EA, sEA, cEA = model.return_EAs()

    model.free()

    return [np.mean(np.abs(EA - EA_in)), np.amax(np.abs(EA - EA_in)),
            np.mean(np.abs(sEA - sEA_in)), np.amax(np.abs(sEA - sEA_in)),
            np.mean(np.abs(cEA - cEA_in)), np.amax(np.abs(cEA - cEA_in))]


if __name__ == "__main__":

    ecc = 1 - 10**(np.linspace(0, -5, 100))
    
    print('#     Ecc       |E - Etrue|        |sinE - sinEtrue|     |cosE - cosEtrue|')
    print('#        ' + 3*'       mean      worst')
    
    for _ecc in ecc:
        errs = calc_accuracy(100000, _ecc)
        print('%9.6f %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e' %
              tuple([_ecc] + errs))
