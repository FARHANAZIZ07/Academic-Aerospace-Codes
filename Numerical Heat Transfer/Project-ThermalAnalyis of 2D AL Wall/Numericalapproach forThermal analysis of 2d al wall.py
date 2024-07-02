# -*- coding: utf-8 -*-
"""
Created on Sat Nov 04 21:39:08 2023

@author: Farhan
"""

import numpy as np
import copy
import seaborn as sns
import matplotlib.pyplot as plt

# Define parameters
k = 200  # W / m K
Lx = 1  # m
Ly = 1  # m
thickness = 0.01  # m
Cp = 903  # J/kg K
rho = 2700  # Kg/m^3
T_air = 200  # K
h_values = [5, 10, 50]  # W/m^2 K
w = 1.1  # over-relaxation w>1
Nx = 32
Ny = 32
Convergence_criteria = 0.00001
delta_x = Lx / (Nx - 2)  # distance between nodes
delta_y = Ly / (Ny - 2)
final_time = 600
delta_t = 1
T_initial = np.full((Ny, Nx, int(final_time/delta_t)), 300)
T_S = 300  # K
T_W = 300  # K
T_E = 300  # K
T_N = 500  # K at t = 0
results = []

for h in h_values:
    T_time = []
    T_time.append(T_initial)
    t_list = [t * delta_t for t in range(0, int(final_time / delta_t))]
    iteration = 0
    T_center_initial = 300
    T_center_values = []
    T_center_values.append(T_center_initial)

    for t in t_list:
        T_guess = T_time[-1]
        converged = False
        while not converged:
            T_new = copy.deepcopy(T_guess)
            for i in range(Nx):
                for j in range(Ny):
                    S_star = h/thickness * (T_new[i,j,t] - T_air)
                    Sp = h/thickness
                    Sc = S_star - Sp * T_new[i,j,t]
                    # Boundary conditions:
                    if i == 0 and 0 <= j <= Nx-1:
                        T_new[i,j,t] = T_N
                    elif 0 <= i < Ny-1 and j == 0:
                        T_new[i,j,t] = T_W
                    elif 0 <= i < Ny-1 and j == Nx-1:
                        T_new[i,j,t] = T_E
                    elif i == Ny-1 and 0 <= j <= Nx-1:
                        T_new[i,j,t] = T_S
                    # Inside nodes
                    elif i == 1 and j == 1:
                        ae = k*delta_y/delta_x
                        aw = k*delta_y/delta_x
                        aN = 2*k*delta_x/delta_y
                        aS = k*delta_x/delta_y
                        ap_0 = rho*Cp*delta_x*delta_y/delta_t
                        b = Sc*delta_x*delta_y + ap_0*T_guess[i,j,t-1]
                        ap = ap_0 + ae + aw + aN + aS - Sp*delta_x*delta_y
                        T_new[i,j,t] = (ae * T_new[i + 1, j,t] + aw * T_new[i - 1, j,t] + aN * T_new[i,j + 1,t] + aS * T_new[i,j - 1,t] + b) / ap
                    elif i == 1 and j == Nx-2:
                        ae = 2*k*delta_y/delta_x
                        aw = k*delta_y/delta_x
                        aN = k*delta_x/delta_y
                        aS = k*delta_x/delta_y
                        ap_0 = rho*Cp*delta_x*delta_y/delta_t
                        b = Sc*delta_x*delta_y + ap_0*T_guess[i,j,t-1]
                        ap = ap_0 + ae + aw + aN + aS - Sp*delta_x*delta_y
                        T_new[i,j,t] = (ae * T_new[i + 1, j,t] + aw * T_new[i - 1, j,t] + aN * T_new[i,j + 1,t] + aS * T_new[i,j - 1,t] + b) / ap
                    elif 1 < i < Ny-2 and j == 1:
                        ae = k*delta_y/delta_x
                        aw = k*delta_y/delta_x
                        aN = k*delta_x/delta_y
                        aS = k*delta_x/delta_y
                        ap_0 = rho*Cp*delta_x*delta_y/delta_t
                        b = Sc*delta_x*delta_y + ap_0*T_guess[i,j,t-1]
                        ap = ap_0 + ae + aw + aN + aS - Sp*delta_x*delta_y
                        T_new[i,j,t] = (ae * T_new[i + 1, j,t] + aw * T_new[i - 1, j,t] + aN * T_new[i,j + 1,t] + aS * T_new[i,j - 1,t] + b) / ap
                    elif i == Ny-2 and j == 1:
                        ae = k*delta_y/delta_x
                        aw = 2*k*delta_y/delta_x
                        aN = k*delta_x/delta_y
                        aS = k*delta_x/delta_y
                        ap_0 = rho*Cp*delta_x*delta_y/delta_t
                        b = Sc*delta_x*delta_y + ap_0*T_guess[i,j,t-1]
                        ap = ap_0 + ae + aw + aN + aS - Sp*delta_x*delta_y
                        T_new[i,j,t] = (ae * T_new[i + 1, j,t] + aw * T_new[i - 1, j,t] + aN * T_new[i,j + 1,t] + aS * T_new[i,j - 1,t] + b) / ap
                    elif 1 < i < Ny-2 and j == Nx-2:
                        ae = 2*k*delta_y/delta_x
                        aw = k*delta_y/delta_x
                        aN = k*delta_x/delta_y
                        aS = k*delta_x/delta_y
                        ap_0 = rho*Cp*delta_x*delta_y/delta_t
                        b = Sc*delta_x*delta_y + ap_0*T_guess[i,j,t-1]
                        ap = ap_0 + ae + aw + aN + aS - Sp*delta_x*delta_y
                        T_new[i,j,t] = (ae * T_new[i + 1, j,t] + aw * T_new[i - 1, j,t] + aN * T_new[i,j + 1,t] + aS * T_new[i,j - 1,t] + b) / ap
                    elif i == Ny-2 and 1 < j < Nx-2:
                        ae = k*delta_y/delta_x
                        aw = k*delta_y/delta_x
                        aN = k*delta_x/delta_y
                        aS = 2*k*delta_x/delta_y
                        ap_0 = rho*Cp*delta_x*delta_y/delta_t
                        b = Sc*delta_x*delta_y + ap_0*T_guess[i,j,t-1]
                        ap = ap_0 + ae + aw + aN + aS - Sp*delta_x*delta_y
                        T_new[i,j,t] = (ae * T_new[i + 1, j,t] + aw * T_new[i - 1, j,t] + aN * T_new[i,j + 1,t] + aS * T_new[i,j - 1,t] + b) / ap
                    elif i == 1 and 1 < j < Nx-2:
                        ae = k*delta_y/delta_x
                        aw = k*delta_y/delta_x
                        aN = 2*k*delta_x/delta_y
                        aS = k*delta_x/delta_y
                        ap_0 = rho*Cp*delta_x*delta_y/delta_t
                        b = Sc*delta_x*delta_y + ap_0*T_guess[i,j,t-1]
                        ap = ap_0 + ae + aw + aN + aS - Sp*delta_x*delta_y
                        T_new[i,j,t] = (ae * T_new[i + 1, j,t] + aw * T_new[i - 1, j,t] + aN * T_new[i,j + 1,t] + aS * T_new[i,j - 1,t] + b) / ap
                    elif i == Ny-2 and j == Nx-2:
                        ae = 2*k*delta_y/delta_x
                        aw = k*delta_y/delta_x
                        aN = k*delta_x/delta_y
                        aS = k*delta_x/delta_y
                        ap_0 = rho*Cp*delta_x*delta_y/delta_t
                        b = Sc*delta_x*delta_y + ap_0*T_guess[i,j,t-1]
                        ap = ap_0 + ae + aw + aN + aS - Sp*delta_x*delta_y
                        T_new[i,j,t] = (ae * T_new[i + 1, j,t] + aw * T_new[i - 1, j,t] + aN * T_new[i,j + 1,t] + aS * T_new[i,j - 1,t] + b) / ap

            # Calculate the residual and check for convergence
            residual = np.abs(T_guess - T_new)
            max_residual = np.max(residual)

            # Update T_new using over-relaxation
            T_new_updated = w * T_new + (1 - w) * T_guess
            T_guess = copy.deepcopy(T_new_updated)  # Update T_guess for the next iteration

            if max_residual <= Convergence_criteria:
                converged = True

            iteration += 1

        Temp_center = T_new_updated[Ny//2, Nx//2, t]
        T_center_values.append(Temp_center)
        T_time.append(T_new_updated)

    # Store the results for the current value of h
    results.append((T_time, T_center_values))

# Plot the initial 2D temperature field for h = 5
sns.set()
plt.figure(figsize=(10, 5))
sns.heatmap(results[0][0][-1][:, :, 10], cmap='plasma', annot=False, cbar=True)
plt.title(f'Temperature at t = 10 seconds for h = {h_values[0]}')
plt.xlabel('Nodes in x-direction')
plt.ylabel('Nodes in y-direction')
plt.show()

sns.set()
plt.figure(figsize=(10, 5))
sns.heatmap(results[0][0][-1][:, :, -1], cmap='plasma', annot=False, cbar=True)
plt.title(f'Temperature at t = 600 seconds for h = {h_values[0]}')
plt.xlabel('Nodes in x-direction')
plt.ylabel('Nodes in y-direction')
plt.show()

# Plot T_center values vs. time for all values of h
plt.figure(figsize=(10, 5))
for i, (T_time, T_center_values) in enumerate(results):
    plt.plot(T_center_values, label=f'T_center (h = {h_values[i]})')

plt.title('Temperature at the Center Node vs. Time')
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.legend()
plt.show()
