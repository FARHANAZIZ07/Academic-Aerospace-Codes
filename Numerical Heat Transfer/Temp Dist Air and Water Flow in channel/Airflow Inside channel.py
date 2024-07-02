# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 19:27:02 2023

@author: Farhan
"""
import numpy as np
import matplotlib.pyplot as plt
# Knowns
L = 0.01 # m
u1 = 0.1 # m/s
u2 = 1 # m/s
velocities = [u1, u2]
Tin = 300 # K
Tout = 400 # K
k = 0.025 # W/mK
rho = 1.226 # kg/m^3
Cp = 1006 # J/kgk
N = 12
x = np. linspace (0, L, N)
# Calculations before the final Equation
gamma = k / Cp
T_profiles_a = [ ]
P_values = []
for u in velocities:
    P = rho *u* L / gamma
    T = Tin + (Tout-Tin) *(np.exp(P * x / L) -1) /(np.exp(P)-1)
    P_values.append(P)
    T_profiles_a.append(T)
    print (f'Peclet number for u = {u} m/s: {P}')
# Plotting the results
for i, u in enumerate(velocities):
    plt.plot(x, T_profiles_a[i], label=f'u = {u} m/s')

plt.xlabel('Distance (m)')
plt.ylabel('Temperature (K)')
plt. legend ()
plt.title(f'Temperature Profile for Different Fluid Velocities for {N} points')
plt.grid(True)
plt.show()


####
delta_x = L / (N - 2) # distance between nodes

# Calculations before Iteration
gamma= k / Cp

# Gauss Seidel Iteration
w = 1 # Coefficient for convergence
Convergence_criteria = 0.0000001 # Convergence criteria

# Create an array to store the temperature profiles for both u1 and u2
T_profiles_b = []

for u in velocities:
    T_guess = np.full(N, 300) # T_Guess
    residual = Convergence_criteria + 1
    T_new = np.zeros(N)
    iteration = 0
    Pe = np.zeros(N)
    Pw =np.zeros(N)
    while residual > Convergence_criteria:
        for i in range(N):
            if i==0:
                ap = 1
                aw = 0
                ae = 0
                b = Tin
                T_new[i] =(ae * T_guess [i + 1] + b) / ap
                Pe[i] = 0
                Pw[i] = 0
            elif i == 1:
                De = gamma/delta_x
                Fe = rho*u
                Dw = 2*gamma / delta_x
                Fw = rho*u
                aw = Dw + max (Fw, 0)
                ae = De + max(-Fe, 0)
                ap = ae + aw + (Fe - Fw)
                b = 0
                T_new[i] = (ae * T_guess [i + 1] + aw * T_guess [i - 1] + b) / ap
                Pe[i] = Fe/De
                Pw[i] = Fw/Dw
            elif 1 < i < N - 2:
                De = gamma/delta_x
                Fe = rho*u
                Dw = gamma / delta_x
                Fw = rho*u
                aw = Dw + max (Fw, 0)
                ae = De + max(-Fe, 0)
                ap = ae + aw + (Fe - Fw)
                b = 0
                T_new[i] = (ae * T_guess [i + 1] + aw * T_guess [i - 1] + b) / ap
                Pe[i] = Fe/De
                Pw[i] = Fw/Dw
            elif i == N - 2:
                De = 2 * gamma/delta_x
                Fe = rho*u
                Dw = gamma / delta_x
                Fw = rho*u
                aw = Dw + max (Fw, 0)
                ae = De + max(-Fe, 0)
                ap = ae + aw + (Fe - Fw)
                b = 0
                T_new[i] = (ae * T_guess [i + 1] + aw * T_guess [i - 1] + b) / ap
                Pe[i] = Fe/De
                Pw[i] = Fw/Dw
            elif i == N - 1:
                ap = 1
                aw = 0
                ae = 0
                b = Tout
                T_new[i]= (aw *T_guess [i-1] + b) /ap
                Pe[i] = 0
                Pw[i] = 0
                    
        residual = max (abs (T_new-T_guess))
        T_new_updated = w * T_new + (1 - w) * T_guess
        T_guess = T_new_updated
        iteration = iteration + 1
        T_profiles_b.append(T_new_updated)
# Plot the temperature profiles for both u1 and u2
x = np. linspace (0, L, N)
plt.plot(x, T_profiles_b[0], label=f'u={u1} m/s')
plt.plot(x, T_profiles_b[1], label=f'u={u2} m/s')
plt.xlabel('x(m)')
plt.ylabel('Temperature (K)')
plt. legend ()
plt.title(f'Temperature Profile for Different Velocities for {N} nodes')
plt.grid(True)
plt.show()



######
fig, ax = plt.subplots()
x = np. linspace (0, L, N)
ax.plot(x, T_profiles_a[0], label=f'u={u} m/s (analytical solution)')
ax.plot(x, T_profiles_a[1], label=f'u={u2} m/s (anlaytical solution)')
ax.plot(x, T_profiles_b[0], label=f'u={u1} m/s (upwind solution)')
ax.plot(x, T_profiles_b[1], label=f'u={u2} m/s (upwind solution)')
ax.set_xlabel('x (m)')
ax.set_ylabel( 'Temperature (K)')
ax. legend()
ax.set_title(f'Combined Temperature Profiles for Different Velocities for {N} Nodes')
ax.grid(True)
plt.show()