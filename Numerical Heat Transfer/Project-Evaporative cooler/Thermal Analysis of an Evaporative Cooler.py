# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 22:35:23 2023

@author: Farhan
"""

import numpy as np
import copy
import matplotlib.pyplot as plt

# Known parameters
k = 26.24*10**(-3)  # W/mK of air
rho = 1.225  # kg/m^3 of air
Cp = 1006  # J/kgK of air
u = 0.5  # m/s
v = 0  # m/s
L=0.02 #m
W=0.003 #m
Lx = L  # m
Ly = W  # m
V_cooler = 0.12*0.12*0.02  # m^3
m_dot_w = 0.02  # kg/s
A_c = L * W #cross section of 2cm by 1cm in m2
X = 0.86 #fraction of the latent energy
rho_w = 997 #kg/m^3 of water
gamma = k / Cp #diffusion coefficient
T_ambient = 294.15 #K

# Source Term calculations for all Latent Energy
V_initial_w = 450 #ml = cm3
V_final_water = 385 #ml = cm3
delta_t_source = 10*60 #sec
V_sponge = V_cooler #m3
delta_V_water = (V_initial_w - V_final_water) / 1000000 #m^3
m_water = delta_V_water * rho_w #kg
Latent_Energy = 2088*1000 #J/kg
E = Latent_Energy * m_water #J
S = E / (V_sponge * delta_t_source) #W/s

# For Gauss-Seidel
w = 1
Nx = 100
Ny = 25
Convergence_criteria = 0.0001
delta_x = Lx / (Nx - 2) # distance between nodes
delta_y = Ly / (Ny - 2)
final_time = 600.0 #10 minutes in seconds
delta_t = 1
T_initial = np.full((Ny, Nx, int(final_time)), 300, dtype = float) #initial guess of temperature profile
t_list = [t*delta_t for t in range(0, int(final_time / delta_t))]
T_time = []
T_time.append(T_initial)
temperature_evolution = []
iteration=0

# For T_bulk and T_out calculations
T_bulk_initial = 357.05
T_bulk = []
dEcv_dt = []
T_time_middle = []
#T_bulk.append(T_bulk_initial)
dEcv_dt=np.full(len(t_list), T_bulk_initial)
T_bulk = np.full(len(t_list), T_bulk_initial)
Q = (1 - X) * S * V_sponge
m_dot_w = 0.02 # kg/s
Cp_w = 4181 # J/kg K
Cv_w = 4138 # J

for t in range(len(t_list)):
    T_guess = T_time[-1]
    converged = False
    T_out= (m_dot_w * Cp_w * T_bulk[t-1] - Q)/ (m_dot_w * Cp_w)
    dEcv_dt[t] = (m_dot_w * Cp_w * T_out - m_dot_w * Cp_w * T_bulk[t-1])
    T_bulk[t] = (dEcv_dt[t] * delta_t)/(m_water*Cv_w) + T_bulk[t-1]
    while not converged:
           T_new = copy.deepcopy(T_guess)
           for i in range(Ny):
               for j in range(Nx):
                    Sp = 0 #because we averaged it
                    Sc = -S/Cp #using negative in here
                    Sc = X*Sc #some energy comes from air some from water so X allows us to know how muc
                    ap_0 = rho*delta_x*delta_y/delta_t #assume rho_p_0 is the same throughout the time
                    #Boundary conditions:
                    if i == 0 and 0 <= j <= Nx-1: #South boundary
                        T_new[i,j] = T_bulk
                    elif 0 < i < Ny-1 and j == 0: #West boundary
                        T_new[i,j] = T_ambient
                    elif 0 < i < Ny-1 and j == Nx-1: #East boundary
                        T_new[i,j] = T_ambient
                    elif i == Ny-1 and 0 <= j <= Nx-1: #North boundary
                        T_new[i,j] = T_bulk
                    #Inside nodes
                    elif i == 1 and j == 1: #South-west corner node
                        De = gamma*delta_y/delta_x
                        Dw = 2*gamma*delta_y/delta_x
                        Dn = gamma*delta_x/delta_y
                        Ds = 2*gamma*delta_x/delta_y
                        Fe = rho*u*delta_y
                        Fw = rho*u*delta_y
                        Fn = rho*v*delta_x
                        Fs = rho*v*delta_x
                        ae = De + max(-Fe, 0)
                        aw = Dw + max(Fw, 0)
                        aS = Ds + max(Fs, 0)
                        aN = Dn + max(-Fn, 0)
                        ap = ae + aw + aN + aS + ap_0 - Sp*delta_x*delta_y
                        b = Sc*delta_x*delta_y + ap_0*T_guess[i,j]
                        T_new[i,j] = (ae * T_new[i,j+1] + aw * T_new[i,j-1] + aN * T_new[i-1,j] + aS * T_new[i+1,j] + b) / ap
                        residual = T_new[i,j] - T_guess[i,j]
                    elif i == 1 and 1 < j < Nx-2: #South central node
                        De = gamma*delta_y/delta_x
                        Dw = gamma*delta_y/delta_x
                        Dn = gamma*delta_x/delta_y
                        Ds = 2*gamma*delta_x/delta_y
                        Fe = rho*u*delta_y
                        Fw = rho*u*delta_y
                        Fn = rho*v*delta_x
                        Fs = rho*v*delta_x
                        ae = De + max(-Fe, 0)
                        aw = Dw + max(Fw, 0)
                        aS = Ds + max(Fs, 0)
                        aN = Dn + max(-Fn, 0)
                        ap = ae + aw + aN + aS + ap_0 - Sp*delta_x*delta_y
                        b = Sc*delta_x*delta_y + ap_0*T_guess[i,j]
                        T_new[i,j] = (ae * T_new[i,j+1] + aw * T_new[i,j-1] + aN * T_new[i-1,j] + aS * T_new[i+1,j] + b) / ap
                        residual = T_new[i,j] - T_guess[i,j]
                    elif i == 1 and j == Nx-2: #South-east corner node
                        De = 2*gamma*delta_y/delta_x
                        Dw = gamma*delta_y/delta_x
                        Dn = gamma*delta_x/delta_y
                        Ds = 2*gamma*delta_x/delta_y
                        Fe = rho*u*delta_y
                        Fw = rho*u*delta_y
                        Fn = rho*v*delta_x
                        Fs = rho*v*delta_x
                        ae = De + max(-Fe, 0)
                        aw = Dw + max(Fw, 0)
                        aS = Ds + max(Fs, 0)
                        aN = Dn + max(-Fn, 0)
                        ap = ae + aw + aN + aS + ap_0 - Sp*delta_x*delta_y
                        b = Sc*delta_x*delta_y + ap_0*T_guess[i,j]
                        T_new[i,j] = (ae * T_new[i,j+1] + aw * T_new[i,j-1] + aN * T_new[i-1,j] + aS * T_new[i+1,j] + b) / ap
                        residual = T_new[i,j] - T_guess[i,j]
                    elif 1 < i < Ny-2 and j == 1: #West central nodes
                        De = gamma*delta_y/delta_x
                        Dw = 2*gamma*delta_y/delta_x
                        Dn = gamma*delta_x/delta_y
                        Ds = gamma*delta_x/delta_y
                        Fe = rho*u*delta_y
                        Fw = rho*u*delta_y
                        Fn = rho*v*delta_x
                        Fs = rho*v*delta_x
                        ae = De + max(-Fe, 0)
                        aw = Dw + max(Fw, 0)
                        aS = Ds + max(Fs, 0)
                        aN = Dn + max(-Fn, 0)
                        ap = ae + aw + aN + aS + ap_0 - Sp*delta_x*delta_y
                        b = Sc*delta_x*delta_y + ap_0*T_guess[i,j]
                        T_new[i,j] = (ae * T_new[i,j+1] + aw * T_new[i,j-1] + aN * T_new[i-1,j] + aS * T_new[i+1,j] + b) / ap
                        residual = T_new[i,j] - T_guess[i,j]
                    elif i == Ny-2 and j == 1: #North-west corner node
                        De = gamma*delta_y/delta_x
                        Dw = 2*gamma*delta_y/delta_x
                        Dn = 2*gamma*delta_x/delta_y
                        Ds = gamma*delta_x/delta_y
                        Fe = rho*u*delta_y
                        Fw = rho*u*delta_y
                        Fn = rho*v*delta_x
                        Fs = rho*v*delta_x
                        ae = De + max(-Fe, 0)
                        aw = Dw + max(Fw, 0)
                        aS = Ds + max(Fs, 0)
                        aN = Dn + max(-Fn, 0)
                        ap = ae + aw + aN + aS + ap_0 - Sp*delta_x*delta_y
                        b = Sc*delta_x*delta_y + ap_0*T_guess[i,j]
                        T_new[i,j] = (ae * T_new[i,j+1] + aw * T_new[i,j-1] + aN * T_new[i-1,j] + aS * T_new[i+1,j] + b) / ap
                        residual = T_new[i,j] - T_guess[i,j]
                    elif i == Ny-2 and 1< j < Nx-2: #North central node
                        De = gamma*delta_y/delta_x
                        Dw = gamma*delta_y/delta_x
                        Dn = 2*gamma*delta_x/delta_y
                        Ds = gamma*delta_x/delta_y
                        Fe = rho*u*delta_y
                        Fw = rho*u*delta_y
                        Fn = rho*v*delta_x
                        Fs = rho*v*delta_x
                        ae = De + max(-Fe, 0)
                        aw = Dw + max(Fw, 0)
                        aS = Ds + max(Fs, 0)
                        aN = Dn + max(-Fn, 0)
                        ap = ae + aw + aN + aS + ap_0 - Sp*delta_x*delta_y
                        b = Sc*delta_x*delta_y + ap_0*T_guess[i,j]
                        T_new[i,j] = (ae * T_new[i,j+1] + aw * T_new[i,j-1] + aN * T_new[i-1,j] + aS * T_new[i+1,j] + b) / ap
                        residual = T_new[i,j] - T_guess[i,j]    
                    elif i == Ny-2 and j == Nx-2: #North-east corner node
                        De = 2*gamma*delta_y/delta_x
                        Dw = gamma*delta_y/delta_x
                        Dn = 2*gamma*delta_x/delta_y
                        Ds = gamma*delta_x/delta_y
                        Fe = rho*u*delta_y
                        Fw = rho*u*delta_y
                        Fn = rho*v*delta_x
                        Fs = rho*v*delta_x
                        ae = De + max(-Fe, 0)
                        aw = Dw + max(Fw, 0)
                        aS = Ds + max(Fs, 0)
                        aN = Dn + max(-Fn, 0)
                        ap = ae + aw + aN + aS + ap_0 - Sp*delta_x*delta_y
                        b = Sc*delta_x*delta_y + ap_0*T_guess[i,j]
                        T_new[i,j] = (ae * T_new[i,j+1] + aw * T_new[i,j-1] + aN * T_new[i-1,j] + aS * T_new[i+1,j] + b) / ap
                        residual = T_new[i,j] - T_guess[i,j]
                    elif 1 < i < Ny-2 and j == Nx-2:  # East central nodes
                        De =2*gamma*delta_y/delta_x
                        Dw = gamma*delta_y/delta_x
                        Dn = gamma*delta_x/delta_y
                        Ds = gamma*delta_x/delta_y
                        Fe = rho*u*delta_y
                        Fw = rho*u*delta_y
                        Fn = rho*v*delta_x
                        Fs = rho*v*delta_x
                        ae = De + max(-Fe, 0)
                        aw = Dw + max(Fw, 0)
                        aS = Ds + max(Fs, 0)
                        aN = Dn + max(-Fn, 0)
                        ap = ae + aw + aN + aS + ap_0 - Sp*delta_x*delta_y
                        b = Sc*delta_x*delta_y + ap_0*T_guess[i,j]
                        T_new[i,j] = (ae * T_new[i,j+1] + aw * T_new[i,j-1] + aN * T_new[i-1,j] + aS * T_new[i+1,j] + b) / ap
                        residual = T_new[i,j] - T_guess[i,j]
                    elif 1 < i < Ny-2 and 1 < j < Nx-2:  # Interior nodes
                        De = gamma*delta_y/delta_x
                        Dw = gamma*delta_y/delta_x
                        Dn = gamma*delta_x/delta_y
                        Ds = gamma*delta_x/delta_y
                        Fe = rho*u*delta_y
                        Fw = rho*u*delta_y
                        Fn = rho*v*delta_x
                        Fs = rho*v*delta_x
                        ae = De + max(-Fe, 0)
                        aw = Dw + max(Fw, 0)
                        aS = Ds + max(Fs, 0)
                        aN = Dn + max(-Fn, 0)
                        ap = ae + aw + aN + aS + ap_0 - Sp*delta_x*delta_y
                        b = Sc*delta_x*delta_y + ap_0*T_guess[i,j]
                        T_new[i,j] = (ae * T_new[i,j+1] + aw * T_new[i,j-1] + aN * T_new[i-1,j] + aS * T_new[i+1,j] + b) / ap
                        residual = T_new[i,j] - T_guess[i,j]
           max_residual = max(residual)
           iteration = iteration+1
           converged = True
           T_guess =copy. deepcopy(T_new)
            
           Temp_middle = T_new[Ny-1, Nx//2] # for Numerical T_out
            
           for i in range (0,Ny):
                 for j in range (0,Nx):
                     if np.max(abs( T_new[i,j]) - abs(T_guess[i,j])) > Convergence_criteria:
                         converged = False
           T_time.append(T_guess)
    T_time_middle = Temp_middle # for numerical T_out   
    temperature_evolution.append(T_new)
    
#To keep same colors for the temperature plots
min_temp = np.min(T_time[-1])
max_temp = np.max(T_time[-1])
x_values = np.linspace(0, Lx, Nx)
y_values = np.linspace(0, Ly, Ny)
X, Y = np.meshgrid(x_values,y_values)

# Colormap oftemperatures at t = 1min
fig = plt.figure(figsize=(10, 5))
temperature_at_1_min = T_new[:,:,60]
plt.pcolormesh(X, Y, temperature_at_1_min, cmap = 'plasma',vmin=min_temp, vmax=max_temp)
cbar = plt.colorbar()
cbar.set_label("Temperature (k)")
plt.xlabel('lenght in x diraction(0 - 0.02 m)')
plt.ylabel('lenght in y diraction(0 - 0.03 m)')
plt.title('Temperature Distribution of Evaporative Cooling channel at t= 1min')
plt.show()

# Colormap oftemperatures at t = 5min
fig = plt.figure(figsize=(10, 5))
temperature_at_5_min = T_new[:,:,300]
plt.pcolormesh(X, Y, temperature_at_5_min, cmap = 'plasma',vmin=min_temp, vmax=max_temp)
cbar = plt.colorbar()
cbar.set_label("Temperature (k)")
plt.xlabel('lenght in x diraction(0 - 0.02 m)')
plt.ylabel('lenght in y diraction(0 - 0.03 m)')
plt.title('Temperature Distribution of Evaporative Cooling channel at t= 5min')
plt.show()

# Colormap oftemperatures at t = 1min
fig = plt.figure(figsize=(10,5))
temperature_at_10_min = T_new[:,:,-1]
plt.pcolormesh(X, Y, temperature_at_10_min, cmap = 'plasma',vmin=min_temp, vmax=max_temp)
cbar = plt.colorbar()
cbar.set_label("Temperature (k)")
plt.xlabel('lenght in x diraction(0 - 0.02 m)')
plt.ylabel('lenght in y diraction(0 - 0.03 m)')
plt.title('Temperature Distribution of Evaporative Cooling channel at t= 10min')
plt.show()



# Plotting
# sns.set()
# plt.figure(figsize=(10, 5))
# sns.heatmap(temperature_evolution[-1][:, :, 60], cmap='plasma', annot=False, cbar=True)
# plt.title('Temperature Distribution at t = 1 min')
# plt.xlabel('Length in x direction (0-0.02 meters)')
# plt.ylabel('Height in y direction (0-0.003 meters)')
# plt.show()

# Comparison of experimental and numerical data
T_bulk_exp = [83.96, 39.7, 30.3, 24.7, 21.2, 19.1, 17.6, 16.6, 15.9, 15.2, 14.8]
T_bulk_exp = [temp + 273.15 for temp in T_bulk_exp]
T_air_out_exp = [21, 30.7, 26.8, 22.7, 20.9, 18.8, 17.8, 16.3, 16.1, 15.4, 15.6]
T_air_out_exp = [temp + 273.15 for temp in T_air_out_exp]
time_exp = [0, 60, 120, 180, 240, 300, 360, 420, 480, 540, 600]

# Plotting experimental data
plt.plot(time_exp, T_bulk_exp, label="Experimental - T_bulk")
plt.plot(time_exp, T_air_out_exp, label="Experimental - T_air_out")
# Plotting numerical data
plt.plot(t_list, T_bulk, label="Numerical T_bulk", linestyle='-', marker='o')
plt.plot(t_list, T_time_middle, label="Numerical T_out", linestyle='-', marker='o')
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.legend()
plt.grid(True)
plt.title("Experimental vs. Numerical Temperature Data")
plt.show()
