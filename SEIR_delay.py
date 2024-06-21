import numpy as np
import matplotlib.pyplot as plt
from my_ddeint import ddeint
import sys, os

POPULATION = 1.4 * (10**9)  # chiny
#POPULATION = 1. # normed to sum 1
C = POPULATION

u_times = []
u_value = []


def u_threshold(i_suma):  # tylko dwa stany - blokada i rozluźnienie
    if i_suma > 100000: 
        return 0.85
    else:
        return 0.0


def u_stopniowe(i_suma):  # stopniowe rozluznianie
    if i_suma > 100000: return 0.85
    if 100000 > i_suma > 50000: return 0.8
    if 50000 >= i_suma > 10000: return 0.5
    if 10000 > i_suma > 100:
        return 0.2
    else:
        return 0.0


def u_brak(i_suma):  # brak reakcji rządu
    return 0.0;


def seir_delay(Y, t, u, tau):
    #stau, itau = Y(t - tau)
    stau, e1tau, e2tau, i1tau, i2tau, i3tau, rtau = Y(t - tau)

    s, e1, e2, i1, i2, i3, r = Y(t)
    gamma = 0.303  #Współczynnik wyzdrowienia (1 / 3.3)
    beta = 2.6 * gamma  # Współczynnik transmisji (R_0 * gamma)
    alpha = 0.196  # Średni okres inkubacji (1 / 5.1)
    mu = 0  # Współczynnik śmiertelności

    #i_suma = i1 + i2 + i3

    i_suma = C*(i1tau + i2tau + i3tau)
    u_times.append(t)    # t might be different than those in ts
    u_value.append(u(i_suma))

    dsdt = -(1 - u(i_suma)) * beta * s * (i1 + i2 + i3) / N # + mu * r0
    de1dt = +(1 - u(i_suma)) * beta * s * (i1 + i2 + i3) / N - 2 * alpha * e1
    de2dt = -2 * alpha * e2                             + 2 * alpha * e1
    di1dt = 2 * alpha * e2 - 3 * gamma * i1 - mu * i1
    di2dt = 3 * gamma * i1 - 3 * gamma * i2 - mu * i2
    di3dt = 3 * gamma * i2 - 3 * gamma * i3 - mu * i3
    drdt = 3 * gamma * i3

    return np.array([dsdt, de1dt, de2dt, di1dt, di2dt, di3dt, drdt])


def initial_history_func(t):
    #return [0.9999, 27504.0, 4538.0, 4538.0, 2.0, 1.0, 2.6]
    #return [9999, 27504.0, 4538.0, 4538.0, 2.0, 1.0, 2.6]  # wartosci z pracy dot chin
    return [0.9999, 0.0001, 0.0, 0.0, 0.0, 0.0, 0.0]  # dla tych wartosci widac najlepsze wyniki


gamma = 0.303  # Współczynnik wyzdrowienia (1 / 3.3)
# N = 1000  # Całkowita populacja
N = 1
R0 = 2.6  # współczynnik reprodukcji
alpha = 0.196  # Średni okres inkubacji (1 / 5.1)

allowed_u = {
    "u_brak": u_brak,
    "u_stopniowe": u_stopniowe,
    "u_threshold": u_threshold,
}

param_tau = 0.0;
if len(sys.argv) > 1:
    param_tau = float(sys.argv[1])

param_u = "";
if len(sys.argv) > 2:
    param_u = sys.argv[2]
if param_u not in allowed_u.keys():
    param_u = "u_brak"

# TODO: replace str+str with os.path.join(...) in file/folder names 
plots_dir = "plots-{}".format(param_u)
os.makedirs(plots_dir, exist_ok=True)
with open(plots_dir + "/.gitignore", "w") as f:
    f.write("*")

print ("param_tau={}, param_u={}".format(param_tau, param_u))

ts = np.linspace(0, 150, 6000) #0, 40, 6000
ys = ddeint(seir_delay, initial_history_func, ts, fargs=(allowed_u[param_u], param_tau,))
#ys = ddeint(seir_delay, initial_history_func, ts, fargs=(u_brak, param_tau,))
#ys = ddeint(seir_delay, initial_history_func, ts, fargs=(u_threshold, param_tau,))
#ys = ddeint(seir_delay, initial_history_func, ts, fargs=(u_stopniowe, param_tau,))

fig, axs = plt.subplots(2, 3, figsize=(11, 6))
fig.suptitle('SEIR Model')

S = C * ys[:, 0]
E = C * (ys[:, 1] + ys[:, 2])
I = C * (ys[:, 3] + ys[:, 4] + ys[:, 5])
R = C * ys[:, 6]

axs[0, 0].plot(ts, S, color='blue')
axs[0, 0].set_title('Podatni (S)')
axs[0, 0].set_xlabel('Czas [dni]')
axs[0, 0].set_ylabel('Populacja')
#axs[0, 0].set_xlim(0, 150)
#axs[0, 0].set_ylim(0, np.max(S) * 1.1)

#"""
axs[0, 1].plot(ts, E, color='orange')
axs[0, 1].set_title('Wystawieni (E)')
axs[0, 1].set_xlabel('Czas [dni]')
axs[0, 1].set_ylabel('Populacja')


axs[1, 0].plot(ts, I, color='red')
axs[1, 0].set_title('Chorzy (I)')
axs[1, 0].set_xlabel('Czas [dni]')
axs[1, 0].set_ylabel('Populacja')

axs[1, 1].plot(ts, R, color='green')
axs[1, 1].set_title('Ozdrowiency (R)')
axs[1, 1].set_xlabel('Czas [dni]')
axs[1, 1].set_ylabel('Populacja')

# Wszystkie
axs[0, 2].plot(ts, S, label='Podatni (S)', color='blue')
axs[0, 2].plot(ts, E, label='Wystawieni (E)', color='orange')
axs[0, 2].plot(ts, I, label='Chorzy (I)', color='red')
axs[0, 2].plot(ts, R, label='Ozdrowiency (R)', color='green')
axs[0, 2].set_title('Model SEIR')
axs[0, 2].set_xlabel('Czas [dni]')
axs[0, 2].set_ylabel('Populacja')
axs[0, 2].legend()
#"""

# value of u w.r.t. t
axs[1, 2].plot(u_times, u_value, color='black')
axs[1, 2].set_title('u(I(t-tau))')
axs[1, 2].set_xlabel('Czas [dni]')
axs[1, 2].set_ylabel('Watosc u')

plt.tight_layout()
#plt.show()
# TODO: replace str+str with os.path.join(...) in file/folder names 
plt.savefig(plots_dir + '/plot-{}.png'.format(param_tau))