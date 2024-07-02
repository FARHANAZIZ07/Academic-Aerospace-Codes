# -*- coding: utf-8 -*-
"""
Created on Fri SEP 22 11:04:18 2023

@author: Farhan
"""

import numpy as np
import matplotlib.pyplot as plt
import copy

#%% Problem #1

convergence_criteria_list = [0.1, 0.001, 0.00001]
T_list = []

for convergence_criteria in convergence_criteria_list:

    # Given Variables and Parameters
    
    k = 5 # W / m K
    L = 1 # m
    q_left = 100 # W/m^2 Heat flux into the left wall
    T_right = 300 # Kelvin. The T of the right hand side
    N = 50 # Number of nodes
    
    # Calculated Variables and Parameters
    
    Deltax = L / (N - 2) # The width of each control volume. 2 Control volumes are infinitely thin
    delx = Deltax
    weight = 1
    
    T_guess = [T_right for item in range(0,N)]
        
    # Enter a while loop that is used to track convergence
    
    converged = False
    
    iterations = 0 # For tracking the number of iterations
    
    while converged == False:
        
        # First start a new list that we will use to track the changes we make to T_guess
        
        T_new = copy.deepcopy(T_guess) # This is important as T_new = T_guess will just create a pointer to T_guess
        iterations = iterations + 1
        x_list = []
    
        # First iterate through each control volume and calculate a new value for T
        # Using the most recent guess value (T_new)
        
        for i in range(0,N):
            
            if i == 0:
                
                ap = k / (delx/2)
                ae = k / (delx/2)
                b = q_left
                x_list.append(0)
                
                T_new[i] = (ae* T_new[i+1] + b) / ap
    
            elif i == 1:
                
                ap = k / (delx/2) + k / (delx)
                aw = k / (delx/2)
                ae = k / (delx)
                b = 0
                x_list.append(delx/2)
                
                T_new[i] = (aw*T_new[i-1] + ae*T_new[i+1] + b) / ap
                
            elif i == N-2:
                
                ap = k / (delx/2) + k / (delx)
                aw = k / (delx)
                ae = k / (delx/2)
                b = 0
                x_list.append(L - delx/2)
                
                T_new[i] = (aw*T_new[i-1] + ae*T_new[i+1] + b) / ap
            
            elif i == N-1:
                
                ap = 1
                b = T_right
                x_list.append(L)

                T_new[i] = b / ap
                
            else:
                
                ap = k / (delx) + k / (delx)
                aw = k / (delx)
                ae = k / (delx)
                b = 0
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
        
    T_list.append(T_new)
    print(iterations)

# Now plot the results

plt.plot(x_list,T_list[0],'k-')
plt.plot(x_list,T_list[1],'r-')
plt.plot(x_list,T_list[2],'b-')
plt.xlabel('Location (m)')
plt.ylabel('Temperature [K]')
plt.legend(['Criteria = 0.1', 'Criteria = 0.001', 'Criteria = 0.00001'])
plt.show()

#%% Problem #2

N_list = [5, 10, 100]
T_list = []
x_list_list = []

# Given Variables and Parameters

k = 5 # W / m K
L = 1 # m
h = 10 # W/m^2 K heat transfer coefficient on the left wall
T_infinity = 400 # Kelvin. The surroundigns temperature on the left wall
T_right = 300 # Kelvin. The T of the right hand side

for N in N_list:
    
    # Calculated Variables and Parameters

    Deltax = L / (N - 2) # The width of each control volume. 2 Control volumes are infinitely thin
    delx = Deltax
    x_list = []
    
    A = np.zeros((N,N))
    b = np.zeros((N,1))
    
    for i in range(0,N):
        
        if i == 0:
            
            A[i][i] = h + k / (delx/2)
            A[i][i+1] = - k / (delx/2)
            b[i] = h*T_infinity
            x_list.append(0)

            
        elif i == 1:
            
            A[i][i] = k / (delx/2) + k / (delx)
            A[i][i-1] = - k / (delx/2)
            A[i][i+1] = - k / (delx)
            b[i] = 0
            x_list.append(delx/2)
            
        elif i == N-2:
            
            A[i][i] = k / (delx/2) + k / (delx)
            A[i][i-1] = - k / (delx)
            A[i][i+1] = - k / (delx/2)
            b[i] = 0
            x_list.append(L - delx/2)
        
        elif i == N-1:
            
            A[i][i] = 1
            b[i] = T_right
            x_list.append(L)
            
        else:
            
            A[i][i] = k / (delx) + k / (delx)
            A[i][i-1] = - k / (delx)
            A[i][i+1] = - k / (delx)
            b[i] = 0
            x_list.append(delx/2 + delx*(i - 1))
            
    # Now solve for the temperature using matrix inversion
    
    T = np.linalg.solve(A,b)
    
    T_list.append(T)
    x_list_list.append(x_list)

# Now plot the results

plt.plot(x_list_list[0],T_list[0],'k-')
plt.plot(x_list_list[1],T_list[1],'r-')
plt.plot(x_list_list[2],T_list[2],'b-')
plt.xlabel('Location (m)')
plt.ylabel('Temperature [K]')
plt.legend(['N = 5', 'N = 10', 'N = 100'])
plt.show()