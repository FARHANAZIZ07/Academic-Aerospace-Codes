# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 17:21:13 2023

@author: Farhan
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import copy

# Given Variables and Parameters

k = 0.6 # W / m K
Cp = 4181 # J / kg K
rho = 1000 # kg/ m3
L = 1000 # m
W = 0.1 # m
T_bottom = 310 # Kelvin
T_top = 310 # kelvin
T_left = 300 # kelvin
T_right = 290 # kevlin
Nx = 200 # The number of nodes in the x direction
Ny = 50 # The number of nodes in the y direction
converged = False

# Velocity Terms

um = 0.1 # m/s
r_o = W/2
r = abs(np.linspace(-r_o, r_o,Ny))
u_list = [2*um*(1-((r_temp/r_o)**2)) for r_temp in r] 

# Calculated Variables and Parameters

Deltax = L / (Nx - 2) # The width of each control volume. 2 Control volumes are infinitely thin
delx = Deltax
Deltay = W / (Ny - 2)
dely = Deltay
position_matrix = []
converged = False
weight = 1.0
convergence_criteria = 0.0001

T_guess = [[T_top for j in range(0,Nx)] for i in range(0,Ny)] 

while converged == False:
    
    T_new = copy.deepcopy(T_guess)
    
    for i in range(0,Ny): # Iterating along y in the outer loop ( i = 0 is the bottom)
        
        for j in range(0,Nx): # Iterating along x in the inner loop ( j = 0 is the left hand side)
        
            if j == 0: # Left side boundary condition
                
                ap = 1
                b = T_left
                T_new[i][j] = b / ap
                
            elif j == Nx - 1: # Right side boundary condition
                
                ap = 1
                b = T_right
                T_new[i][j] = b / ap
                
            elif i == 0: # Bottom boundary condition
                
                ap = 1
                b = T_bottom
                T_new[i][j] = b / ap
                
            elif i == Ny - 1: # Top side boundary condition
                
                ap = 1
                b = T_top
                T_new[i][j] = b / ap
                
            elif i == 1 and j == 1: # Bottom left inner control volume
                
                u = u_list[i]
                v = 0
                
                # Diffusion Terms
            
                De = (k*Deltay)/(Cp*delx)
                Dw = (k*Deltay)/(Cp*delx/2)
                Dn = (k*Deltax)/(Cp*dely)
                Ds = (k*Deltax)/(Cp*dely/2)
                
                # Advection Terms
            
                Fe = rho*u*Deltay
                Fw = rho*u*Deltay
                Fn = rho*v*Deltax
                Fs = rho*v*Deltax
                
                aw = Dw + max([Fw,0])
                ae = De + max([-Fe,0])
                an = Dn + max([-Fn,0])
                as1 = Ds + max([Fs,0])
                ap = (aw + ae + an + as1)
                b = 0
                T_new[i][j] = (ae*T_new[i][j+1] + aw*T_new[i][j-1] + an*T_new[i+1][j] + as1*T_new[i-1][j] + b) / ap
                        
            elif i == 1 and j == Nx - 2: # Bottom right inner control volume
                
                u = u_list[i]
                v = 0
                
                # Diffusion Terms
            
                De = (k*Deltay)/(Cp*delx/2)
                Dw = (k*Deltay)/(Cp*delx)
                Dn = (k*Deltax)/(Cp*dely)
                Ds = (k*Deltax)/(Cp*dely/2)
                
                # Advection Terms
            
                Fe = rho*u*Deltay
                Fw = rho*u*Deltay
                Fn = rho*v*Deltax
                Fs = rho*v*Deltax
                
                aw = Dw + max([Fw,0])
                ae = De + max([-Fe,0])
                an = Dn + max([-Fn,0])
                as1 = Ds + max([Fs,0])
                ap = (aw + ae + an + as1)
                b = 0
                T_new[i][j] = (ae*T_new[i][j+1] + aw*T_new[i][j-1] + an*T_new[i+1][j] + as1*T_new[i-1][j] + b) / ap
    
            elif i == Ny - 2 and j == Nx - 2: # Top right inner control volume
                
                u = u_list[i]
                v = 0
                
                # Diffusion Terms
            
                De = (k*Deltay)/(Cp*delx/2)
                Dw = (k*Deltay)/(Cp*delx)
                Dn = (k*Deltax)/(Cp*dely/2)
                Ds = (k*Deltax)/(Cp*dely)
                
                # Advection Terms
            
                Fe = rho*u*Deltay
                Fw = rho*u*Deltay
                Fn = rho*v*Deltax
                Fs = rho*v*Deltax
                
                aw = Dw + max([Fw,0])
                ae = De + max([-Fe,0])
                an = Dn + max([-Fn,0])
                as1 = Ds + max([Fs,0])
                ap = aw + ae + an + as1
                b = 0
                T_new[i][j] = (ae*T_new[i][j+1] + aw*T_new[i][j-1] + an*T_new[i+1][j] + as1*T_new[i-1][j] + b) / ap
            
            elif i == Ny - 2 and j == 1: # Top left inner control volume
                
                u = u_list[i]
                v = 0
                
                # Diffusion Terms
            
                De = (k*Deltay)/(Cp*delx)
                Dw = (k*Deltay)/(Cp*delx/2)
                Dn = (k*Deltax)/(Cp*dely/2)
                Ds = (k*Deltax)/(Cp*dely)
                
                # Advection Terms
            
                Fe = rho*u*Deltay
                Fw = rho*u*Deltay
                Fn = rho*v*Deltax
                Fs = rho*v*Deltax
                
                aw = Dw + max([Fw,0])
                ae = De + max([-Fe,0])
                an = Dn + max([-Fn,0])
                as1 = Ds + max([Fs,0])
                ap = aw + ae + an + as1
                b = 0
                T_new[i][j] = (ae*T_new[i][j+1] + aw*T_new[i][j-1] + an*T_new[i+1][j] + as1*T_new[i-1][j] + b) / ap
    
            elif j == 1: # Any other CV in the inner left column
                
                u = u_list[i]
                v = 0
                
                # Diffusion Terms
            
                De = (k*Deltay)/(Cp*delx)
                Dw = (k*Deltay)/(Cp*delx/2)
                Dn = (k*Deltax)/(Cp*dely)
                Ds = (k*Deltax)/(Cp*dely)
                
                # Advection Terms
            
                Fe = rho*u*Deltay
                Fw = rho*u*Deltay
                Fn = rho*v*Deltax
                Fs = rho*v*Deltax
                
                aw = Dw + max([Fw,0])
                ae = De + max([-Fe,0])
                an = Dn + max([-Fn,0])
                as1 = Ds + max([Fs,0])
                ap = aw + ae + an + as1
                b = 0
                T_new[i][j] = (ae*T_new[i][j+1] + aw*T_new[i][j-1] + an*T_new[i+1][j] + as1*T_new[i-1][j] + b) / ap
            
            elif j == Nx - 2: # Any other CV in the inner right column
                
                u = u_list[i]
                v = 0
                
                # Diffusion Terms
            
                De = (k*Deltay)/(Cp*delx/2)
                Dw = (k*Deltay)/(Cp*delx)
                Dn = (k*Deltax)/(Cp*dely)
                Ds = (k*Deltax)/(Cp*dely)
                
                # Advection Terms
            
                Fe = rho*u*Deltay
                Fw = rho*u*Deltay
                Fn = rho*v*Deltax
                Fs = rho*v*Deltax
                
                aw = Dw + max([Fw,0])
                ae = De + max([-Fe,0])
                an = Dn + max([-Fn,0])
                as1 = Ds + max([Fs,0])
                ap = aw + ae + an + as1
                b = 0
                T_new[i][j] = (ae*T_new[i][j+1] + aw*T_new[i][j-1] + an*T_new[i+1][j] + as1*T_new[i-1][j] + b) / ap
            
            elif i == 1: # Any other CV in the inner bottom row
                
                u = u_list[i]
                v = 0
                
                # Diffusion Terms
            
                De = (k*Deltay)/(Cp*delx)
                Dw = (k*Deltay)/(Cp*delx)
                Dn = (k*Deltax)/(Cp*dely)
                Ds = (k*Deltax)/(Cp*dely/2)
                
                # Advection Terms
            
                Fe = rho*u*Deltay
                Fw = rho*u*Deltay
                Fn = rho*v*Deltax
                Fs = rho*v*Deltax
                
                aw = Dw + max([Fw,0])
                ae = De + max([-Fe,0])
                an = Dn + max([-Fn,0])
                as1 = Ds + max([Fs,0])
                ap = aw + ae + an + as1
                b = 0
                T_new[i][j] = (ae*T_new[i][j+1] + aw*T_new[i][j-1] + an*T_new[i+1][j] + as1*T_new[i-1][j] + b) / ap
            
            elif i == Ny - 2: # Any other CV in the inner top row
                
                u = u_list[i]
                v = 0
                
                # Diffusion Terms
            
                De = (k*Deltay)/(Cp*delx)
                Dw = (k*Deltay)/(Cp*delx)
                Dn = (k*Deltax)/(Cp*dely/2)
                Ds = (k*Deltax)/(Cp*dely)
                
                # Advection Terms
            
                Fe = rho*u*Deltay
                Fw = rho*u*Deltay
                Fn = rho*v*Deltax
                Fs = rho*v*Deltax
                
                aw = Dw + max([Fw,0])
                ae = De + max([-Fe,0])
                an = Dn + max([-Fn,0])
                as1 = Ds + max([Fs,0])
                ap = aw + ae + an + as1
                b = 0
                T_new[i][j] = (ae*T_new[i][j+1] + aw*T_new[i][j-1] + an*T_new[i+1][j] + as1*T_new[i-1][j] + b) / ap
            
            else:
                
                u = u_list[i]
                v = 0
                
                # Diffusion Terms
            
                De = (k*Deltay)/(Cp*delx)
                Dw = (k*Deltay)/(Cp*delx)
                Dn = (k*Deltax)/(Cp*dely)
                Ds = (k*Deltax)/(Cp*dely)
                
                # Advection Terms
            
                Fe = rho*u*Deltay
                Fw = rho*u*Deltay
                Fn = rho*v*Deltax
                Fs = rho*v*Deltax
                
                aw = Dw + max([Fw,0])
                ae = De + max([-Fe,0])
                an = Dn + max([-Fn,0])
                as1 = Ds + max([Fs,0])
                ap = aw + ae + an + as1
                b = 0
                T_new[i][j] = (ae*T_new[i][j+1] + aw*T_new[i][j-1] + an*T_new[i+1][j] + as1*T_new[i-1][j] + b) / ap
                
#            print('ae: ' + str(ae))
#            print('aw: ' + str(aw))

    converged = True
    convergence_list = [[0 for j in range(0,Nx)] for i in range(0,Ny)]
    
    for i in range(0,Ny):
        
        for j in range(0,Nx):
            
            convergence_list[i][j] = (abs(abs(T_new[i][j]) - abs(T_guess[i][j])))
            
            if abs(abs(T_new[i][j]) - abs(T_guess[i][j])) > convergence_criteria:
                
                converged = False
                
                #print('repeating')
    
    T_guess = [[weight*T_new[i][j] + (1 - weight)*T_guess[i][j] for j in range(0,Nx)] for i in range(0,Ny)] 

    #sns.heatmap(convergence_list)
    #plt.show()

# Now plot the results

#plt.figure(figsize = (L*30,W*30))
sns.heatmap(T_guess, square=True)
plt.show()