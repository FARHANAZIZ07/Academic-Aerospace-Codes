# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 21:18:33 2024

@author: Farhan
"""

import random
import numpy as np
import matplotlib.pyplot as plt

# Define constants
w = 1  # m # Defining width of surfaces
L = 1  # m # Constants for analytical solution
w_i = 1  # m # Constants for analytical solution
w_j = 1  # m # Constants for analytical solution
L_values = np.arange(0.5, 5.5, 0.5)  # Defining length between surfaces
N_rays_values = [1000, 10000, 100000]  # Defining number of rays emitted, needs to be large

# Initialize an empty list to store Fij values analytical
Fij_values = []

# Loop to do it over different N_rays values
for N_rays in N_rays_values:
    view_factor_values = []  # Initialize an empty list to store view factors for numerical
    Fij_values_temp = []  # Initialize an empty list to store Fij values analytical

    # Loop to do it over different L values
    for L in L_values:
        # Analytical calculations from equation given in part 1
        W_i = w_i / L
        W_j = w_j / L
        Fij = (((W_i + W_j)**2 + 4)**(1/2) - ((W_j - W_i)**2 + 4)**(1/2)) / (2 * W_i)
        Fij_values_temp.append(Fij)

        # Numerical calculations:
        # Creating Dictionary for surface 1 and 2 where we will store the rays to keep track
        Surf1 = {'x0_ray': [], 'y0_ray': [], 'xhat_ray': [], 'yhat_ray': [],
                 'emission_points': [], 'emission_vector': [], 'tr': [], 'ts': [],
                 'Impacted_surfaces': [], 'emission': N_rays, 'Rt_values': []}

        Surf2 = {'x0_ray': [], 'y0_ray': [], 'xhat_ray': [], 'yhat_ray': [],
                 'emission_points': [], 'emission_vector': [], 'tr': [], 'ts': [],
                 'Impacted_surfaces': [], 'emission': N_rays, 'Rt_values': []}

        # Setting parametric vectors from drawing in notes, and Eq (1)
        Surf1['x0_s'] = 0  # Setting x0 for surface 1 for parametric vector
        Surf1['x1_s'] = w  # Setting x1 for surface 1 for parametric vector
        Surf1['y0_s'] = 0  # Setting y0 for surface 1 for parametric vector
        Surf1['y1_s'] = 0  # Setting y1 for surface 1 for parametric vector
        Surf2['x0_s'] = w  # Setting x0 for surface 2 for parametric vector
        Surf2['x1_s'] = 0  # Setting x1 for surface 2 for parametric vector
        Surf2['y0_s'] = L  # Setting y0 for surface 2 for parametric vector
        Surf2['y1_s'] = L  # Setting y1 for surface 2 for parametric vector

        # Calculating xhat and yhat for parametric vector for each surface from Eq (2)
        Surf1['xhat_s'] = Surf1['x1_s'] - Surf1['x0_s']  # xhat Surface 1
        Surf2['xhat_s'] = Surf2['x1_s'] - Surf2['x0_s']  # xhat Surface 2
        Surf1['yhat_s'] = Surf1['y1_s'] - Surf1['y0_s']  # yhat Surface 1
        Surf2['yhat_s'] = Surf2['y1_s'] - Surf2['y0_s']  # yhat Surface 2

        # Calculating normal vectors to the surface vectors from Eq (4)
        Surf1['xhat90'] = -Surf1['yhat_s']  # xhat of 90 degrees for Surface 1
        Surf1['yhat90'] = Surf1['xhat_s']  # yhat of 90 degrees for Surface 1
        Surf2['xhat90'] = -Surf2['yhat_s']  # xhat of 90 degrees for Surface 2
        Surf2['yhat90'] = Surf2['xhat_s']  # yhat of 90 degrees for Surface 2

        Surf1['Impacted_surfaces'] = []

        # Calculating Emission Points and Directions
        for i in range(N_rays):
            # Gather Vector Data from Directory needed, in this problem only emitting surface 1
            x0_s = Surf1['x0_s']
            xhat_s = Surf1['xhat_s']
            y0_s = Surf1['y0_s']
            yhat_s = Surf1['yhat_s']
            xhat90 = Surf1['xhat90']
            yhat90 = Surf1['yhat90']

            # Generating Random number between 0 and 1 to generate a random spot for each ray
            Rt = random.random()
            ts = Rt
            Surf1['Rt_values'].append(Rt)  # Saving the random number to keep track

            # Calculate x0 and y0 (Location of ray emitted) and save to Directory
            x0_ray = x0_s + xhat_s * ts  # Calculate x0_ray from Eq (3)
            y0_ray = y0_s + yhat_s * ts  # Calculate y0_ray from Eq (3)
            Surf1['x0_ray'].append(x0_ray)  # Save to Surface 1 Directory
            Surf1['y0_ray'].append(y0_ray)  # Save to Surface 1 Directory
            Surf2['x0_ray'].append(x0_ray)  # Save to Surface 2 Directory - Not needed just in case
            Surf2['y0_ray'].append(y0_ray)  # Save to Surface 2 Directory - Not needed just in case
            Surf1['emission_points'].append([x0_ray, y0_ray])  # Save to Surface 1 Directory as points
            Surf2['emission_points'].append([x0_ray, y0_ray])  # Save to Surface 2 Directory as points

            # Generate random theta and phi to find random emission angle
            R_theta = random.random()  # Random number between 0 and 1
            theta = np.arcsin(R_theta)  # Calculate theta from random number Eq (5)
            R_phi = random.random()  # Random number between 0 and 1

            # Above theta will always be positive so we also want to account
            # for negative so by generating random number we set some theta negative
            if R_phi < 0.5:
                theta = -theta
            else:
                theta = theta

            # Calculate ray pointing vector (xhat, yhat) and Save to Directory
            xhat_ray = np.cos(theta) * xhat90 - np.sin(theta) * yhat90  # From Eq (6)
            yhat_ray = np.sin(theta) * xhat90 + np.cos(theta) * yhat90  # From Eq (6)
            Surf1['xhat_ray'].append(xhat_ray)  # Save to Surface 1 Directory
            Surf1['yhat_ray'].append(yhat_ray)  # Save to Surface 1 Directory
            Surf2['xhat_ray'].append(xhat_ray)  # Save to Surface 2 Directory - Not needed just in case
            Surf2['yhat_ray'].append(yhat_ray)  # Save to Surface 2 Directory - Not needed just in case
            Surf1['emission_vector'].append([xhat_ray, yhat_ray])  # Save to Surface 1 Directory as points
            Surf2['emission_vector'].append([xhat_ray, yhat_ray])  # Save to Surface 2 Directory as points

        # Calculating Intersections between rays emitted from surface 1 and surface 2
        for i in range(N_rays):
            # Gather Vector Data from Directory needed, so receiving rays is surface 2
            x0_s = Surf2['x0_s']
            y0_s = Surf2['y0_s']
            xhat_s = Surf2['xhat_s']
            yhat_s = Surf2['yhat_s']

            # Gather Vector Data from Directory needed, so emitting rays is surface 1
            x0_ray = Surf1['x0_ray'][i]
            y0_ray = Surf1['y0_ray'][i]
            xhat_ray = Surf1['xhat_ray'][i]
            yhat_ray = Surf1['yhat_ray'][i]

            # Calculating tr from Eq (7) by setting vector eq of surf2 and eq of rays from surf1 equals
            tr = (xhat_s * y0_s - y0_ray * xhat_s + yhat_s * x0_s - xhat_s * x0_ray) / (yhat_ray * xhat_s - yhat_s * xhat_ray)
            ts = (x0_ray + xhat_ray * tr - x0_s) / xhat_s

            # Saving tr and ts into Directory to keep track and verify it is correct
            Surf2['tr'].append(tr)
            Surf2['ts'].append(ts)

            # Determine if the ray impacted or not. If ts is between 0 and 1, it did impact
            # If it impacts, then find the surface with smallest tr, in this case not needed
            # If ts is less than 0 or bigger than 1 then no impact -> Inf
            if ts < 0 or ts > 1:
                Surf2['Impacted_surfaces'].append('Inf')  # Saving non-impact rays as Inf
            elif 0 <= ts <= 1:
                Surf2['Impacted_surfaces'].append(tr)  # Saving impact rays as tr

        # Calculating Number of Impacted rays by counting entries in Impacted_Surfaces that are not Inf
        # meaning the amount of rays that hit a surface
        num_impacted_rays = sum(1 for surface in Surf2['Impacted_surfaces'] if surface != 'Inf')

        # Calculating Total number of rays emitted (should be N_rays as well)
        total_rays = len(Surf2['Impacted_surfaces'])

        # Calculate the view factor F from Eq (10) and save into list to plot
        view_factor = num_impacted_rays / total_rays
        view_factor_values.append(view_factor)

    # Store Fij values analytical
    Fij_values.append(Fij_values_temp)

    # Plotting analytical and numerical solutions
    if N_rays == 1000:
        plt.plot(L_values, Fij_values_temp, marker='o', label='Analytical Solution')
        plt.plot(L_values, view_factor_values, marker='o', label=f'Number of Rays = {N_rays}')

# Show the plot with title and axis titles and grid
plt.legend()
plt.title('View Factor for Numerical and Analytical solutions as a Function of L for Different Number of Rays')
plt.xlabel('L (meters)')
plt.ylabel('View Factors (F_1_2)')
plt.grid(True)
plt.show()
