# Lista de Exercícios de CPC777 - Modelagem e Simulação de Reservatórios
# Problema 2: Pressão no fundo no poço (Pwd)

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
gamma = 1.78108


# Definição de funções adimensionalizadas

def p(pd):
    import math
    return Pi - pd * qostd * Bo * muo / (2 * math.pi * k * h)


def td(t):
    return k * t / (phi * muo * ct * rw ** 2)


# Criando o vetor de tempos para saída de 0.01 s até 10 anos.

tempos = np.linspace(0.01, 10 *365 * 24 * 3600, 200)

# Criando as variáveis admensionais:

tempos_adm = td(tempos)

# Criando os diversos fatores Skin

skin = np.array([-5, -2, 0, 2, 5])

# Avaliar a Equação 1.1 para os tempos e raios escolhidos

ns = len(skin)
nt = len(tempos_adm)

# Criando um array vazio de P e Pd

P = np.zeros((ns, nt))
Pd = np.zeros((ns, nt))


for i in range(0, ns):
    for j in range(0, nt):
        Pd[i][j] = 0.5 * np.log(4 / gamma * tempos_adm[j]) + skin[i]
        P[i][j] = p(Pd[i][j])

# Convertendo P de Pa para psi:
P_psi = P * 1.45038e-4

# Criando saídas Gráficas
plt.semilogx(tempos[:], P_psi[0][:], label="Skin = -5")
plt.semilogx(tempos[:], P_psi[1][:], label="Skin = -2")
plt.semilogx(tempos[:], P_psi[2][:], label="Skin = 0")
plt.semilogx(tempos[:], P_psi[3][:], label="Skin = 2")
plt.semilogx(tempos[:], P_psi[4][:], label="Skin = 5")
plt.xlabel('tempo (s)')
plt.ylabel('Pressão do fundo do poço (psi)')
plt.legend()
plt.show()