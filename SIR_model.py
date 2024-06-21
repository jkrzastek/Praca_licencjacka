import sys
sys.path.append(".")

import matplotlib.pyplot as plt
import numpy as np
from diffsolve.mysolvers import rksolver


def sir_derivs(t, y, beta, alpha):  # SIR without vital dynamics
    s, i, r = y
    dsdt = -beta * s * i
    didt = beta * s * i - alpha * i
    drdt = alpha * i
    return np.array([dsdt, didt, drdt])


# zmieniłem troche parametry, bo nie było widać żadnych zmian.
i0 = 0.001
beta = 3
alpha = 1.5
s0 = 1 - i0
e0 = 0
r0 = 0
y0 = np.array([s0, i0, r0])
dt = 0.1 # dt z doświadczenia - lepiej, żeby było <1, zazwyczaj używa się czegoś w okolicach 0.1, 0.01,
# a jesli wyniki dziwnie się zachowują, funkcje nie wyglądają na gładkie, to mozna spróbowac zmniejszyć dalej.
t_final = 100

# tak się nie uda, bo rksolver oczekuje funkjci derivs z dwoma parametrami, a dostaje z czterema
#result_t, result_y = rksolver(sir_derivs, y0, dt, t_final)

# mozna temu zaradzic, tworząc funkcję, która nadpisuje część parametrów, na te dostepne
# w tej chwili w przestrzeni nazw, najlepiej do tego nadaje się konstrukcja lambda
# która jest w zasadzie tym samym co def nazwa_funkcji, ale może być wpisana w jednej linijce:
sir_with_parameters_set = lambda t, y: sir_derivs(t, y, beta, alpha)
# i potem uzywamy juz tej funkcji
print(sir_with_parameters_set(0, y0))
result_t, result_y = rksolver(sir_with_parameters_set, y0, dt, t_final) 

plt.plot(result_t, result_y[:, 0], label='Susceptible')
plt.plot(result_t, result_y[:, 1], label='Infectious')
plt.plot(result_t, result_y[:, 2], label='Recovered')
plt.xlabel('Time')
plt.xlim(0,22)
plt.ylabel('Population')
plt.title('SIR model solution')
plt.legend()
plt.show()
