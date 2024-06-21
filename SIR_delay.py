import numpy as np
import matplotlib.pyplot as plt

def rksolver_with_delay(derivs, y, dt, t_final, tau1, tau2):
    time = 0
    Nsteps = round(t_final / dt)
    t = np.zeros((Nsteps, 1))
    data = np.zeros((Nsteps, y.shape[0]))
    t[0] = time
    data[0, :] = y

    for i in range(1, Nsteps):
        k1 = dt * derivs(time, y, tau1, tau2)
        k2 = dt * derivs(time + dt/2, y + k1/2, tau1, tau2)
        k3 = dt * derivs(time + dt/2, y + k2/2, tau1, tau2)
        k4 = dt * derivs(time + dt, y + k3, tau1, tau2)
        y = y + k1/6 + k2/3 + k3/3 + k4/6
        time = time + dt
        t[i] = time
        data[i, :] = y
    return t, data

def sir_derivs_with_delay(t, y, tau1, tau2):
    s, i, r = y
    dsdt = -beta * s * i
    didt = beta * s * i - alpha * i
    drdt = alpha * i
    if t > tau1:
        didt -= alpha * (i - np.interp(t - tau1, t_history, init_I_history))
    if t > tau2:
        drdt -= gamma * (r - np.interp(t - tau2, t_history, init_R_history))
    return np.array([dsdt, didt, drdt])


i0 = 0.001
beta = 3
alpha = 0.5
gamma = 0.001
s0 = 1 - i0
r0 = 0
y0 = np.array([s0, i0, r0])
dt = 0.1
t_final = 100

# Opóźnienia
tau1 = 0 # opóźnienie S->I w dniach od zarażenia do objawów, incubation period
tau2 = 0 # opóźnienie I->R w dniach od zarażenia do wyzdrowienia

# Dane historii inicjalizacyjnych
t_history = np.linspace(0, t_final, int(t_final / dt) + 1)
init_I_history = np.zeros_like(t_history)
init_R_history = np.zeros_like(t_history)

# Inicjalizacja historii dla zainfekowanych (I) i ozdrowiałych (R)
init_I_history[0] = i0
init_R_history[0] = r0

# Wywołanie solvera z opóźnieniem
sir_with_delay_parameters_set = lambda t, y, tau1, tau2: sir_derivs_with_delay(t, y, tau1, tau2)
result_t, result_y = rksolver_with_delay(sir_with_delay_parameters_set, y0, dt, t_final, tau1, tau2)


plt.plot(result_t, result_y[:, 0], label='Susceptible')
plt.plot(result_t, result_y[:, 1], label='Infectious')
plt.plot(result_t, result_y[:, 2], label='Recovered')
plt.xlim(0, 30)
plt.xlabel('Time [days]')
plt.ylabel('Population [%]')
plt.title('SIR model with delay')
plt.legend()
plt.show()
