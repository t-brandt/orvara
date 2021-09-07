from orvara import orbit
import numpy as np
import time


def random_rv_data_params_model(Nsamples):
    jit, mpri, msec, sau = 0.5, 1, 0.1, 10
    esino, ecoso, inc, asc, lam = 0.5, 0.5, 1, 1, 1
    data = orbit.Data(27321, 'None', 'None', 'None')
    data.custom_epochs(list(np.linspace(data.refep, data.refep+1000, Nsamples)))
    model = orbit.Model(data)
    theta = np.array([jit, mpri, msec, sau, esino, ecoso, inc, asc, lam])
    params = orbit.Params(theta, 0, 1, data.nInst, 1)
    orbit.calc_EA_RPP(data, params, model)
    return data, params, model


def time_calc_rv(Nsamples, loops):
    data, params, model = random_rv_data_params_model(Nsamples=int(Nsamples))

    times = []
    for i in range(int(loops)):
        t = time.time_ns()
        orbit.calc_RV(data, params, model)
        times.append(time.time_ns() - t)
    return np.mean(times)/Nsamples, '+-', np.std(times)/np.sqrt(loops), 'nano seconds per eval'


if __name__ == "__main__":
    print(time_calc_rv(1E2, 1E5))