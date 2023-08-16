# Lista de Exercícios de CPC777 - Modelagem e Simulação de Reservatórios
# Problema 1: fluxo radial unidimensional em regime transiente

# Importação de módulos

import matplotlib.pyplot as plt
import numpy as np

# Definição de parâmetros:

phi = 0.2       # fração
ct = 2.18e-9    # Pa
k = 148e-15     # m^2
h = 3.048       # m
qostd = 518e-5  # m^3 s^-1
Bo = 1.5        # m^3 mstd^-3
muo = 3.3e-4    # Pa s
Pi = 15.17e6    # Pa
rw = 0.0889     # m
re = 609.6      # m


# Definição de funções adimensionalizadas

def p(pd):
    import math
    return Pi - pd * qostd * Bo * muo / (2 * math.pi * k * h)


def td(t):
    return k * t / (phi * muo * ct * rw ** 2)


def rd(r):
    return r / rw


# Definição da Função Exponencial Integral, pois a definição da expi() do Scipy é diferente daquela que conhecemos.
# Mais detalhes podem ser consultados em:
# https://transientgroundwaterflow.readthedocs.io/en/latest/TransientFlowToAWell.html
# Fonte da função expint():
# http://www.analyticelements.org/mw/index.php?title=Exponential_Integral_(a.k.a._Theis%27_Well_Function)

def expint(x):
    import math
    w = 0
    if x < 15:
        n = 1
        termo = -x
        w = -0.5772156649 - math.log(x) - termo
        while abs(termo) > 1e-9:
            termo = -x * n / (n+1) ** 2 * termo
            w = w - termo
            n = n + 1
    return w


# Criando o vetor de tempos para saída em 1min, 1 hora, 1 dias, 30 dias, 1 ano e 10 anos.


tempos = np.array([1 * 60, 1 * 3600, 24 * 3600, 30 * 24 * 3600, 365 * 24 * 3600, 10 * 365 * 24 * 3600])

# Criando o vetor de raios variando entre rw e re de 200 pontos

raios = np.linspace(rw, re, 200)

# Criando as variáveis admensionais:

tempos_adm = td(tempos)
raios_adm = rd(raios)

# Avaliar a Equação 1.1 para os tempos e raios escolhidos

nt = len(tempos_adm)
nr = len(raios_adm)

# Criando um array vazio de P e Pd

P = np.zeros((nt, nr))
Pd = np.zeros((nt, nr))


for i in range(0, nt):
    for j in range(0, nr):
        Pd[i][j] = 0.5 * expint(raios_adm[j] ** 2 / (4 * tempos_adm[i]))
        P[i][j] = p(Pd[i][j])

# Convertendo P de Pa para psi:
P_psi = P * 1.45038e-4

# Criando saídas Gráficas
plt.plot(raios[:], P_psi[0][:], label="Tempo = 1 minuto")
plt.plot(raios[:], P_psi[1][:], label="Tempo = 1 hora")
plt.plot(raios[:], P_psi[2][:], label="Tempo = 1 dia")
plt.plot(raios[:], P_psi[3][:], label="Tempo = 30 dias")
plt.plot(raios[:], P_psi[4][:], label="Tempo = 1 ano")
plt.plot(raios[:], P_psi[5][:], label="Tempo = 10 anos")
plt.xlabel('Raio (m)')
plt.ylabel('Pressão (psi)')
plt.legend()
plt.show()
