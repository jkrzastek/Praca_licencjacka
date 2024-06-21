import numpy as np
import matplotlib.pyplot as plt
from my_ddeint import ddeint


def sis(Y, t):
    s, i = Y(t)
    stau, itau = Y(t - 1)
    beta = 1.8  # Wspolczynnik transmisji
    gamma = 0.6  # wspolczynnik wyzdrowien
    dsdt = -beta * s * i + gamma * i
    didt = beta * s * i - gamma * i
    return np.array([dsdt, didt])


def beta(i): #stopniowe rozluznianie
    if i > 0.25: return 0.85
    if 0.25 > i > 0.20: return 1.2
    if 0.20 >= i > 0.15: return 1.4
    if 0.15 > i > 0.1: return 1.6
    else: return 1.85

i_values = np.linspace(0, 0.3, 100)  # Zakres wartości i od 0 do 0.3
beta_values = [beta(i) for i in i_values]  # Wyliczenie wartości beta dla każdego i

plt.plot(i_values, beta_values, label='Funkcja beta(i)')
plt.xlabel('Wartość i')
plt.ylabel('Wartość beta')
plt.title('Wykres funkcji beta(i)')
plt.legend()
#plt.grid(True)
plt.show()

def sis_delay1(Y, t):
    s, i = Y(t)
    stau, itau = Y(t - 1)
    #beta = 1.85  # Wspolczynnik transmisji
    gamma = 0.6  # wspolczynnik wyzdrowien
    # dla R = beta/gamma > 1 choroba rozwija sie

    dsdt = -beta(itau) * stau * itau + gamma * itau
    didt = beta(itau) * stau * itau - gamma * itau
    return np.array([dsdt, didt])


def sis_delay(Y, t, tau):
    stau, itau = Y(t - tau)
    s1, i1 = Y(t - 2)
    s0, i0 = Y(t)
    gamma = 0.6

    dsdt = -beta(itau) * s0 * i0 + gamma * i1
    didt = beta(itau) * s0 * i0 - gamma * i1
    return np.array([dsdt, didt])

def initial_history_func(t):
    return [.999, 0.001]


#itau_values = np.array([0.01, 0.5, 1, 1.25, 1.5, 2])
itau_values = np.array([0.01, 0.01, 0.03, 0.04, 0.05, 0.1])

fig, axs = plt.subplots(2, 3, figsize=(10, 6))
fig.suptitle('SIS dla roznych opoznien itau')


for idx, itau in enumerate(itau_values):
    print("Value of itau:", itau)
    ts = np.linspace(0, 30, 6000)
    ys = ddeint(sis_delay, initial_history_func, ts, fargs=(itau,))

    row = idx // 3
    col = idx % 3

    axs[row, col].plot(ts, ys[:, 0], color='blue', label='Susceptible')
    axs[row, col].plot(ts, ys[:, 1], color='orange', label='Infectious')
    axs[row, col].set_title(f'itau = {itau}')
    axs[row, col].set_xlabel('Czas [dni]')
    axs[row, col].set_ylabel('Populacja')
    axs[row, col].legend()

plt.tight_layout()
plt.show()
