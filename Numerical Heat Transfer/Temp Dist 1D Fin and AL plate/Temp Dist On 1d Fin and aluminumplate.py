# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 11:26:28 2024

@author: Farhan
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import copy

#%% Problem #1 - Fin Problem

# Given Variables and Parameters

k = 401 # W / m K
L = .02 # m
D = 0.003 # Diameter
h = 10 # W / m2 K
T_base = 400 # W/m^2 Heat flux into the left wall
emissivity_list = [1]
T_infinity = 273 # Kelvin. The T of the right hand side
N = 100 # Number of nodes
sigma = 5.67*10**(-8)

# Calculated Variables and Parameters

Deltax = L / (N - 2) # The width of each control volume. 2 Control volumes are infinitely thin
delx = Deltax
weight = 0.5

T_guess = [T_infinity for item in range(0,N)]
    
# Enter a while loop that is used to track convergence

converged = False
convergence_criteria = 0.000001

iterations = 0 # For tracking the number of iterations

for emissivity in emissivity_list:
    
    while converged == False:
        
        # First start a new list that we will use to track the changes we make to T_guess
        
        T_new = copy.deepcopy(T_guess) # This is important as T_new = T_guess will just create a pointer to T_guess
        iterations = iterations + 1
        x_list = []
    
        # First iterate through each control volume and calculate a new value for T
        # Using the most recent guess value (T_new)
        
        for i in range(0,N):
            
            if i == 0:
                
                ap = 1
                b = T_base
                x_list.append(0)
                
                T_new[i] = (b) / ap
    
            elif i == 1:
                
                ap = k / (delx/2) + k / (delx) + (4*h/D)*delx
                aw = k / (delx/2)
                ae = k / (delx)
                b = ((4/D) * (h*T_infinity - emissivity*sigma*(T_new[i]**4 - T_infinity**4)))*delx/2
                x_list.append(delx/2)
                
                T_new[i] = (aw*T_new[i-1] + ae*T_new[i+1] + b) / ap
                
            elif i == N-2:
                
                ap = k / (delx/2) + k / (delx) + (4*h/D)*delx
                aw = k / (delx)
                ae = k / (delx/2)
                b = ((4/D) * (h*T_infinity - emissivity*sigma*(T_new[i]**4 - T_infinity**4)))*delx/2
                x_list.append(L - delx/2)
                
                T_new[i] = (aw*T_new[i-1] + ae*T_new[i+1] + b) / ap
            
            elif i == N-1:
                
                ap = k / (delx/2)
                aw = k / (delx/2)
                b = 0
                x_list.append(L)
    
                T_new[i] = (b + aw*T_new[i-1]) / ap
                
            else:
                
                ap = k / (delx) + k / (delx) + (4*h/D)*delx
                aw = k / (delx)
                ae = k / (delx)
                b = ((4/D) * (h*T_infinity - emissivity*sigma*(T_new[i]**4 - T_infinity**4)))*delx
                x_list.append(delx/2 + delx*(i - 1))
                
                T_new[i] = (aw*T_new[i-1] + ae*T_new[i+1] + b) / ap
                
                
        # Now iterate through the list again and check for convergence. We start by
        # setting converged to be true and then will see if any nodes violate
        # the convergence criteria which will cause us to set it to false
    
        converged = True
        
        for i in range(0,N):
            
            if abs(T_new[i] - T_guess[i]) > convergence_criteria:
                
                converged = False
                
        # Now that we have checked for convergence, we will either repeat the loop
        # again (if converged is still false) or exit the loop. Either way, we will
        # reset T_guess to be our most recent value
        
        T_guess = [weight*T_new[i] + (1 - weight)*T_guess[i] for i in range(0,N)]
        
# Now plot the results
        
    if emissivity == 0:

        plt.plot(x_list,T_new,'k-')
        
    elif emissivity == 1:
        
        plt.plot(x_list,T_new,'r-')
        
plt.xlabel('Location (m)')
plt.ylabel('Temperature [K]')
plt.legend(['Emissivity = 1'])
plt.show()

#%% Problem #2

# Given Variables and Parameters

k = 300 # W / m K
L = 1 # m
W = 1 # m
T_bottom = 400 # Kelvin
T_top = 300 # kelvin
T_left = 300 # kelvin
T_right = 300 # kevlin
Nx = 50 # The number of nodes in the x direction
Ny = 50 # The number of nodes in the y direction
converged = False
    
# Calculated Variables and Parameters

Deltax = L / (Nx - 2) # The width of each control volume. 2 Control volumes are infinitely thin
delx = Deltax
Deltay = W / (Ny - 2)
dely = Deltay
position_matrix = []
converged = False
weight = 1.5
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
                
                aw = k / (delx/2)
                ae = k / (delx)
                an = k / (dely)
                as1 = k / (dely/2)
                ap = (aw + ae + an + as1)
                b = 0
                T_new[i][j] = (ae*T_new[i][j+1] + aw*T_new[i][j-1] + an*T_new[i+1][j] + as1*T_new[i-1][j] + b) / ap
                        
            elif i == 1 and j == Nx - 2: # Bottom right inner control volume
                
                aw = k / (delx)
                ae = k / (delx/2)
                an = k / (dely)
                as1 = k / (dely/2)
                ap = (aw + ae + an + as1)
                b = 0
                T_new[i][j] = (ae*T_new[i][j+1] + aw*T_new[i][j-1] + an*T_new[i+1][j] + as1*T_new[i-1][j] + b) / ap
    
            elif i == Ny - 2 and j == Nx - 2: # Top right inner control volume
                
                aw = k / (delx)
                ae = k / (delx/2)
                an = k / (dely/2)
                as1 = k / (dely)
                ap = aw + ae + an + as1
                b = 0
                T_new[i][j] = (ae*T_new[i][j+1] + aw*T_new[i][j-1] + an*T_new[i+1][j] + as1*T_new[i-1][j] + b) / ap
            
            elif i == Ny - 2 and j == 1: # Top left inner control volume
                
                aw = k / (delx/2)
                ae = k / (delx)
                an = k / (dely/2)
                as1 = k / (dely)
                ap = aw + ae + an + as1
                b = 0
                T_new[i][j] = (ae*T_new[i][j+1] + aw*T_new[i][j-1] + an*T_new[i+1][j] + as1*T_new[i-1][j] + b) / ap
    
            elif j == 1: # Any other CV in the inner left column
                
                aw = k / (delx/2)
                ae = k / (delx)
                an = k / (dely)
                as1 = k / (dely)
                ap = aw + ae + an + as1
                b = 0
                T_new[i][j] = (ae*T_new[i][j+1] + aw*T_new[i][j-1] + an*T_new[i+1][j] + as1*T_new[i-1][j] + b) / ap
            
            elif j == Nx - 2: # Any other CV in the inner right column
                
                aw = k / (delx)
                ae = k / (delx/2)
                an = k / (dely)
                as1 = k / (dely)
                ap = aw + ae + an + as1
                b = 0
                T_new[i][j] = (ae*T_new[i][j+1] + aw*T_new[i][j-1] + an*T_new[i+1][j] + as1*T_new[i-1][j] + b) / ap
            
            elif i == 1: # Any other CV in the inner bottom row
                
                aw = k / (delx)
                ae = k / (delx)
                an = k / (dely)
                as1 = k / (dely/2)
                ap = aw + ae + an + as1
                b = 0
                T_new[i][j] = (ae*T_new[i][j+1] + aw*T_new[i][j-1] + an*T_new[i+1][j] + as1*T_new[i-1][j] + b) / ap
            
            elif i == Ny - 2: # Any other CV in the inner top row
                
                aw = k / (delx)
                ae = k / (delx)
                an = k / (dely/2)
                as1 = k / (dely)
                ap = aw + ae + an + as1
                b = 0
                T_new[i][j] = (ae*T_new[i][j+1] + aw*T_new[i][j-1] + an*T_new[i+1][j] + as1*T_new[i-1][j] + b) / ap
            
            else:
                
                aw = k / (delx)
                ae = k / (delx)
                an = k / (dely)
                as1 = k / (dely)
                ap = aw + ae + an + as1
                b = 0
                T_new[i][j] = (ae*T_new[i][j+1] + aw*T_new[i][j-1] + an*T_new[i+1][j] + as1*T_new[i-1][j] + b) / ap

    converged = True
    convergence_list = [[0 for j in range(0,Nx)] for i in range(0,Ny)]
    
    for i in range(0,Ny):
        
        for j in range(0,Nx):
            
            convergence_list[i][j] = (abs(abs(T_new[i][j]) - abs(T_guess[i][j])))
            
            if abs(abs(T_new[i][j]) - abs(T_guess[i][j])) > convergence_criteria:
                
                converged = False
    
    T_guess = [[weight*T_new[i][j] + (1 - weight)*T_guess[i][j] for j in range(0,Nx)] for i in range(0,Ny)] 

    #sns.heatmap(convergence_list)
    #plt.show()

# Now plot the results

sns.heatmap(T_guess)
plt.show()
    
