# Lista de Exercícios de CPC777 - Modelagem e Simulação de Reservatórios
# Problema 3: Solução do Regime Permanente para o Esquema de Injeção Five-Spot

# Importação de módulos

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

# Definição dos parâmetros

phi = 0.2           # fração
k = 148e-15         # m^2
mu_o = 3.3e-4       # Pa*s
ct = 2.18e-9        # Pa^-1
q_prod = 518e-5     # m^3/s
Bo = 1.5            # m^3/m^3 padrão
Pi = 15.17e6        # Pa
rw = 0.0889         # m
d = 609.6           # m
h = 3.048           # m
n = 50              # número de pontos na malha
m = 50              # número de termos no somatório

# Alocar com zeros para agilizar a execução

x = y = np.linspace(0, d, n)
P = np.zeros((n, n), dtype=float)

# Criando a constante do termo fonte

c = (16 * q_prod * Bo * mu_o / (np.pi ** 2 * k * h))


# Criando o termo do somatório

def somar(rmax, smax, xx, yy, aa, bb):
    a = np.zeros((rmax, smax), dtype=float)
    for r in range(0, rmax):
        for s in range(0, smax):
            a[r, s] = np.cos((2 * r + 1) * np.pi * xx / aa) * np.cos((2 * s + 1) * np.pi * yy / bb) / (
                        (2 * r + 1) ** 2 + (2 * s + 1) ** 2)
    b = np.sum(a)
    return b


# Criando o loop aninhado

for i in range(0, n):
    for j in range(0, n):
        P[i, j] = Pi - c * somar(m, m, x[i], y[j], d, d)

# Convertendo P de Pa para Kgf/cm2

P_psi = P * 1.45038e-4

# Criando saída gráfica

X, Y = np.meshgrid(x, y)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
Z = P_psi.reshape(X.shape)
ax.plot_surface(X, Y, Z, cmap=cm.jet)
# inicial = np.ones((200, 200), dtype=float) * Pi * 1.45038e-4
# ZZ = inicial.reshape(X.shape)
# ax.plot_surface(X, Y, ZZ, alpha=0.2, cmap=cm.winter)

ax.set_xlabel('Direção x (m)')
ax.set_ylabel('Direção y (m)')
ax.set_zlabel('Pressão (psi)')

plt.show()
