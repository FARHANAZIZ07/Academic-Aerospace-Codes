# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 21:11:13 2024

@author: Farhan
"""

import random
import pandas as pd
import numpy as np
# Define constants
L_surface = 0.1 #m
N_surf = 10
w = L_surface/N_surf # m # Defining width of surfaces
V_groove_angle = [30, 90, 120] #degrees
V_groove_angle_half_values = np.deg2rad([angle / 2 for angle in V_groove_angle])
N_rays_values = [500000] # Defining number of rays emitted, needs to be large so after testing I decided to to 500,000
N_surfaces = 1 #Number of surfaces emitting rays
# Initialize an empty list to store Fij values analytical
Fij_values = []

rows=[]

# Loop to do it over different N_rays values
for N_rays in N_rays_values:
    view_factor_values_1_2 = [] # Initialize an empty list to store view factors for numerical
    view_factor_values_1_3 = [] # Initialize an empty list to store view factors for numerical
    view_factor_values_1_4 = [] # Initialize an empty list to store view factors for numerical
    view_factor_values_1_5 = [] # Initialize an empty list to store view factors for numerical
    view_factor_values_1_6 = [] # Initialize an empty list to store view factors for numerical
    view_factor_values_1_7 = [] # Initialize an empty list to store view factors for numerical
    view_factor_values_1_8 = [] # Initialize an empty list to store view factors for numerical
    view_factor_values_1_9 = [] # Initialize an empty list to store view factors for numerical
    view_factor_values_1_10 = [] # Initialize an empty list to store view factors for numerical
    view_factor_values_1_11 = [] # Initialize an empty list to store view factors for numerical
    
    # Loop to do it over different angles values
    for V_groove_angle_half in V_groove_angle_half_values:
        #Numerical calculations:
        # Creating Dictionary for surface 1 (emitting) and 2 - 11 (recieving) where we will store the rays to keep track
        Surf1 = {'x0,ray': [], 'y0,ray': [], 'xhat,ray': [], 'yhat,ray': [],
        'emission_points': [], 'emission_vector': [], 'tr': [], 'ts': [],
        'Impacted_surfaces': [], 'emission': N_rays, 'Rt_values': []}
        Surf2 = {'x0,ray': [], 'y0,ray': [], 'xhat,ray': [], 'yhat,ray': [],
        'emission_points': [], 'emission_vector': [], 'tr': [], 'ts': [],
        'Impacted_surfaces': [], 'emission': N_rays, 'Rt_values': []}
        Surf3 = {'x0,ray': [], 'y0,ray': [], 'xhat,ray': [], 'yhat,ray': [],
        'emission_points': [], 'emission_vector': [], 'tr': [], 'ts': [],
        'Impacted_surfaces': [], 'emission': N_rays, 'Rt_values': []}
        Surf4 = {'x0,ray': [], 'y0,ray': [], 'xhat,ray': [], 'yhat,ray': [],
        'emission_points': [], 'emission_vector': [], 'tr': [], 'ts': [],
        'Impacted_surfaces': [], 'emission': N_rays, 'Rt_values': []}
        Surf5 = {'x0,ray': [], 'y0,ray': [], 'xhat,ray': [], 'yhat,ray': [],
        'emission_points': [], 'emission_vector': [], 'tr': [], 'ts': [],
        'Impacted_surfaces': [], 'emission': N_rays, 'Rt_values': []}
        Surf6 = {'x0,ray': [], 'y0,ray': [], 'xhat,ray': [], 'yhat,ray': [],
        'emission_points': [], 'emission_vector': [], 'tr': [], 'ts': [],
        'Impacted_surfaces': [], 'emission': N_rays, 'Rt_values': []}
        Surf7 = {'x0,ray': [], 'y0,ray': [], 'xhat,ray': [], 'yhat,ray': [],
        'emission_points': [], 'emission_vector': [], 'tr': [], 'ts': [],
        'Impacted_surfaces': [], 'emission': N_rays, 'Rt_values': []}
        Surf8 = {'x0,ray': [], 'y0,ray': [], 'xhat,ray': [], 'yhat,ray': [],
        'emission_points': [], 'emission_vector': [], 'tr': [], 'ts': [],
        'Impacted_surfaces': [], 'emission': N_rays, 'Rt_values': []}
        Surf9 = {'x0,ray': [], 'y0,ray': [], 'xhat,ray': [], 'yhat,ray': [],
        'emission_points': [], 'emission_vector': [], 'tr': [], 'ts': [],
        'Impacted_surfaces': [], 'emission': N_rays, 'Rt_values': []}
        Surf10 = {'x0,ray': [], 'y0,ray': [], 'xhat,ray': [], 'yhat,ray': [],
        'emission_points': [], 'emission_vector': [], 'tr': [], 'ts': [],
        'Impacted_surfaces': [], 'emission': N_rays, 'Rt_values': []}
        Surf11 = {'x0,ray': [], 'y0,ray': [], 'xhat,ray': [], 'yhat,ray': [],
        'emission_points': [], 'emission_vector': [], 'tr': [], 'ts': [],
        'Impacted_surfaces': [], 'emission': N_rays, 'Rt_values': []}
        # Setting parametric vectors from drawing in notes, and Eq (1)
        Surf1['x0,s'] = - 1*w*np.sin(V_groove_angle_half) #Setting x0 for surface 1 for parametric vector
        Surf1['x1,s'] = - 0*w*np.sin(V_groove_angle_half) #Setting x1 for surface 1 for parametric vector
        Surf1['y0,s'] = 1*w*np.cos(V_groove_angle_half) #Setting y0 for surface 1 for parametric vector
        Surf1['y1,s'] = 0*w*np.cos(V_groove_angle_half) #Setting y1 for surface 1 for parametric vector
        Surf2['x0,s'] = 0 #Setting x0 for surface 2 for parametric vector
        Surf2['x1,s'] = w*np.sin(V_groove_angle_half) #Setting x1 for surface 2 for parametric vector
        Surf2['y0,s'] = 0 #Setting y0 for surface 2 for parametric vector
        Surf2['y1,s'] = w*np.cos(V_groove_angle_half) #Setting y1 for surface 2 for parametric vector
        Surf3['x0,s'] = Surf2['x1,s'] #Setting x0 for surface 3 for parametric vector
        Surf3['x1,s'] = 2*w*np.sin(V_groove_angle_half) #Setting x1 for surface 3 for parametric vector
        Surf3['y0,s'] = Surf2['y1,s'] #Setting y0 for surface 3 for parametric vector
        Surf3['y1,s'] = 2*w*np.cos(V_groove_angle_half) #Setting y1 for surface 3 for parametric vector
        Surf4['x0,s'] = Surf3['x1,s'] #Setting x0 for surface 4 for parametric vector
        Surf4['x1,s'] = 3*w*np.sin(V_groove_angle_half) #Setting x1 for surface 4 for parametric vector
        Surf4['y0,s'] = Surf3['y1,s'] #Setting y0 for surface 4 for parametric vector
        Surf4['y1,s'] = 3*w*np.cos(V_groove_angle_half) #Setting y1 for surface 4 for parametric vector
        Surf5['x0,s'] = Surf4['x1,s'] #Setting x0 for surface 5 for parametric vector
        Surf5['x1,s'] = 4*w*np.sin(V_groove_angle_half) #Setting x1 for surface 5 for parametric vector
        Surf5['y0,s'] = Surf4['y1,s'] #Setting y0 for surface 5 for parametric vector
        Surf5['y1,s'] = 4*w*np.cos(V_groove_angle_half) #Setting y1 for surface 5 for parametric vector
        Surf6['x0,s'] = Surf5['x1,s'] #Setting x0 for surface 6 for parametric vector
        Surf6['x1,s'] = 5*w*np.sin(V_groove_angle_half) #Setting x1 for surface 6 for parametric vector
        Surf6['y0,s'] = Surf5['y1,s'] #Setting y0 for surface 6 for parametric vector
        Surf6['y1,s'] = 5*w*np.cos(V_groove_angle_half) #Setting y1 for surface 6 for parametric vector
        Surf7['x0,s'] = Surf6['x1,s'] #Setting x0 for surface 7 for parametric vector
        Surf7['x1,s'] = 6*w*np.sin(V_groove_angle_half) #Setting x1 for surface 7 for parametric vector
        Surf7['y0,s'] = Surf6['y1,s'] #Setting y0 for surface 7 for parametric vector
        Surf7['y1,s'] = 6*w*np.cos(V_groove_angle_half) #Setting y1 for surface 7 for parametric vector
        Surf8['x0,s'] = Surf7['x1,s'] #Setting x0 for surface 8 for parametric vector
        Surf8['x1,s'] = 7*w*np.sin(V_groove_angle_half) #Setting x1 for surface 8 for parametric vector
        Surf8['y0,s'] = Surf7['y1,s'] #Setting y0 for surface 8 for parametric vector
        Surf8['y1,s'] = 7*w*np.cos(V_groove_angle_half) #Setting y1 for surface 8 for parametric vector
        Surf9['x0,s'] = Surf8['x1,s'] #Setting x0 for surface 9 for parametric vector
        Surf9['x1,s'] = 8*w*np.sin(V_groove_angle_half) #Setting x1 for surface 9 for parametric vector
        Surf9['y0,s'] = Surf8['y1,s'] #Setting y0 for surface 9 for parametric vector
        Surf9['y1,s'] = 8*w*np.cos(V_groove_angle_half) #Setting y1 for surface 9 for parametric vector
        Surf10['x0,s'] = Surf9['x1,s'] #Setting x0 for surface 10 for parametric vector
        Surf10['x1,s'] = 9*w*np.sin(V_groove_angle_half) #Setting x1 for surface 10 for parametric vector
        Surf10['y0,s'] = Surf9['y1,s'] #Setting y0 for surface 10 for parametric vector
        Surf10['y1,s'] = 9*w*np.cos(V_groove_angle_half) #Setting y1 for surface 10 for parametric vector
        Surf11['x0,s'] = Surf10['x1,s'] #Setting x0 for surface 11 for parametric vector
        Surf11['x1,s'] = 10*w*np.sin(V_groove_angle_half) #Setting x1 for surface 11 for parametric vector
        Surf11['y0,s'] = Surf10['y1,s'] #Setting y0 for surface 11 for parametric vector
        Surf11['y1,s'] = 10*w*np.cos(V_groove_angle_half) #Setting y1 for surface 11 for parametric vector
        # Calculating xhat and yhat for parametric vector for each surface from Eq (2)
        Surf1['xhat,s'] = Surf1['x1,s'] - Surf1['x0,s'] #xhat Surface 1
        Surf1['yhat,s'] = Surf1['y1,s'] - Surf1['y0,s'] #yhat Surface 1
        Surf2['xhat,s'] = Surf2['x1,s'] - Surf2['x0,s'] #xhat Surface 2
        Surf2['yhat,s'] = Surf2['y1,s'] - Surf2['y0,s'] #yhat Surface 2
        Surf3['xhat,s'] = Surf3['x1,s'] - Surf3['x0,s'] #xhat Surface 3
        Surf3['yhat,s'] = Surf3['y1,s'] - Surf3['y0,s'] #yhat Surface 3
        Surf4['xhat,s'] = Surf4['x1,s'] - Surf4['x0,s'] #xhat Surface 4
        Surf4['yhat,s'] = Surf4['y1,s'] - Surf4['y0,s'] #yhat Surface 4
        Surf5['xhat,s'] = Surf5['x1,s'] - Surf5['x0,s'] #xhat Surface 5
        Surf5['yhat,s'] = Surf5['y1,s'] - Surf5['y0,s'] #yhat Surface 5
        Surf6['xhat,s'] = Surf6['x1,s'] - Surf6['x0,s'] #xhat Surface 6
        Surf6['yhat,s'] = Surf6['y1,s'] - Surf6['y0,s'] #yhat Surface 6
        Surf7['xhat,s'] = Surf7['x1,s'] - Surf7['x0,s'] #xhat Surface 7
        Surf7['yhat,s'] = Surf7['y1,s'] - Surf7['y0,s'] #yhat Surface 7
        Surf8['xhat,s'] = Surf8['x1,s'] - Surf8['x0,s'] #xhat Surface 8
        Surf8['yhat,s'] = Surf8['y1,s'] - Surf8['y0,s'] #yhat Surface 8
        Surf9['xhat,s'] = Surf9['x1,s'] - Surf9['x0,s'] #xhat Surface 9
        Surf9['yhat,s'] = Surf9['y1,s'] - Surf9['y0,s'] #yhat Surface 9
        Surf10['xhat,s'] = Surf10['x1,s'] - Surf10['x0,s'] #xhat Surface 10
        Surf10['yhat,s'] = Surf10['y1,s'] - Surf10['y0,s'] #yhat Surface 10
        Surf11['xhat,s'] = Surf11['x1,s'] - Surf11['x0,s'] #xhat Surface 11
        Surf11['yhat,s'] = Surf11['y1,s'] - Surf11['y0,s'] #yhat Surface 11
        # Calculating normal vectors to the surface vectors from Eq (4)
        Surf1['xhat90'] = -Surf1['yhat,s'] #xhat of 90degrees for Surface 1
        Surf1['yhat90'] = Surf1['xhat,s'] #yhat of 90degrees for Surface 1
        Surf2['xhat90'] = -Surf2['yhat,s'] #xhat of 90degrees for Surface 2
        Surf2['yhat90'] = Surf2['xhat,s'] #yhat of 90degrees for Surface 2
        Surf3['xhat90'] = -Surf3['yhat,s'] #xhat of 90degrees for Surface 3
        Surf3['yhat90'] = Surf3['xhat,s'] #yhat of 90degrees for Surface 3
        Surf4['xhat90'] = -Surf4['yhat,s'] #xhat of 90degrees for Surface 4
        Surf4['yhat90'] = Surf4['xhat,s'] #yhat of 90degrees for Surface 4
        Surf5['xhat90'] = -Surf5['yhat,s'] #xhat of 90degrees for Surface 5
        Surf5['yhat90'] = Surf5['xhat,s'] #yhat of 90degrees for Surface 5
        Surf6['xhat90'] = -Surf6['yhat,s'] #xhat of 90degrees for Surface 6
        Surf6['yhat90'] = Surf6['xhat,s'] #yhat of 90degrees for Surface 6
        Surf7['xhat90'] = -Surf7['yhat,s'] #xhat of 90degrees for Surface 7
        Surf7['yhat90'] = Surf7['xhat,s'] #yhat of 90degrees for Surface 7
        Surf8['xhat90'] = -Surf8['yhat,s'] #xhat of 90degrees for Surface 8
        Surf8['yhat90'] = Surf8['xhat,s'] #yhat of 90degrees for Surface 8
        Surf9['xhat90'] = -Surf9['yhat,s'] #xhat of 90degrees for Surface 9
        Surf9['yhat90'] = Surf9['xhat,s'] #yhat of 90degrees for Surface 9
        Surf10['xhat90'] = -Surf10['yhat,s'] #xhat of 90degrees for Surface 10
        Surf10['yhat90'] = Surf10['xhat,s'] #yhat of 90degrees for Surface 10
        Surf11['xhat90'] = -Surf11['yhat,s'] #xhat of 90degrees for Surface 11
        Surf11['yhat90'] = Surf11['xhat,s'] #yhat of 90degrees for Surface 11
        Surf1['Impacted_surfaces'] = [] #Setting list to store impacted surfaces
        # Calculating Emission Points and Directions
        # Do for every surface emitting rays - so 1 in this case
        
        for i in range(0, 1):
            for j in range(0, N_rays):
                # Gather Vector Data from Directory needed, in this problem only emitting surface 1
                x0_s = Surf1['x0,s']
                xhat_s = Surf1['xhat,s']
                y0_s = Surf1['y0,s']
                yhat_s = Surf1['yhat,s']
                xhat90 = Surf1['xhat90']
                yhat90 = Surf1['yhat90']
                # Generating Random number between 0 and 1 to generate a random spot for each ray
                Rt = random.random()
                ts = Rt
                Surf1['Rt_values'].append(Rt) # Saving the random number to keep track
                # Calculate x0 and y0 (Location of ray emitted) and save to Directory
                x0_ray = x0_s + xhat_s * ts # Calculate x0_ray from Eq (3)
                y0_ray = y0_s + yhat_s * ts # Calculate y0_ray from Eq (3)
                Surf1['x0,ray'].append(x0_ray) #Save to Surface 1 Directory
                Surf1['y0,ray'].append(y0_ray) #Save to Surface 1 Directory
                Surf1['emission_points'].append([x0_ray, y0_ray]) #Save to Surface 1 Directory as points
                # Generate random theta and phi to find random emission angle
                R_theta = random.random() # Random number between 0 and 1
                theta = np.arcsin((R_theta)) # Calculate theta from random number Eq (5)
                R_phi = random.random() # Random phi between 0 and 1
                # Above theta will always be positive so we also want to account
                # for negative so by generating random number we set some theta negative
                if R_phi < 0.5:
                    theta = -theta
                else:
                    theta = theta
                # Calculate ray pointing vector (xhat yhat) and Save to Directory
                xhat_ray = np.cos(theta)*xhat90 - np.sin(theta)*yhat90 # From eq (6)
                yhat_ray = np.sin(theta)*xhat90 + np.cos(theta)*yhat90 # From eq (6)
                Surf1['xhat,ray'].append(xhat_ray) #Save to Surface 1 Directory
                Surf1['yhat,ray'].append(yhat_ray) #Save to Surface 1 Directory
                Surf1['emission_vector'].append([xhat_ray, yhat_ray]) # Save to Surface 1 Directory as points
                
        # Calculating Intersections between rays emitted from surface 1 and surface 2 - 11
        # Do for every surface recieving rays
        for i in range(0, N_surfaces):
            for j in range(0, N_rays):
                for k in range(0, N_surfaces):
                    # Gather Vector Data from Directory needed, so emitting rays is surface 1
                    x0_ray = Surf1['x0,ray'][j]
                    y0_ray = Surf1['y0,ray'][j]
                    xhat_ray = Surf1['xhat,ray'][j]
                    yhat_ray = Surf1['yhat,ray'][j]
                    # Gather Vector Data from Directory needed, so recieveing rays is surface 2
                    x0_s2 = Surf2['x0,s']
                    y0_s2 = Surf2['y0,s']
                    xhat_s2 = Surf2['xhat,s']
                    yhat_s2 = Surf2['yhat,s']
                    # Calculating tr from Eq (7) by setting vector eq of surf2 and eq of rays from surf1 equals
                    tr2 = (xhat_s2 * y0_s2 - y0_ray * xhat_s2 + yhat_s2 * x0_ray - yhat_s2 * x0_s2) / (
                    yhat_ray * xhat_s2 - yhat_s2 * xhat_ray)
                    # Calculating ts from Eq (8) by setting vector eq of surf2 and eq of rays from surf1 equals
                    ts2 = (x0_ray + xhat_ray * tr2 - x0_s2) / xhat_s2
                    # Saving tr and ts into Directory to keep track and verify it is correct
                    Surf2['tr'].append(tr2)
                    Surf2['ts'].append(ts2)
                    # Determine if the ray impacted or not. If ts is between 0 and 1, it did impact
                    # If it impact, then find surface with smallest tr
                    # If ts is less than 0 or bigger than 1 then no impact -> Inf
                    # Comes from Eq (9)
                    if ts2 < 0 or ts2 > 1:
                        Surf2['Impacted_surfaces'].append('Inf') # Saving non impact rays as Inf
                    elif 0 <= ts2 <=1:
                        Surf2['Impacted_surfaces'].append(tr2) # Saving impact rays as tr
                    # Gather Vector Data from Directory needed, so receiving rays is surface 3
                    x0_s3 = Surf3['x0,s']
                    y0_s3 = Surf3['y0,s']
                    xhat_s3 = Surf3['xhat,s']
                    yhat_s3 = Surf3['yhat,s']
                    # Calculating tr from Eq (7) by setting vector eq of surf3 and eq of rays from surf1 equals
                    tr3 = (xhat_s3 * y0_s3 - y0_ray * xhat_s3 + yhat_s3 * x0_ray - yhat_s3 * x0_s3) / (
                    yhat_ray * xhat_s3 - yhat_s3 * xhat_ray)
                    # Calculating ts from Eq (8) by setting vector eq of surf3 and eq of rays from surf1 equals
                    ts3 = (x0_ray + xhat_ray * tr3 - x0_s3) / xhat_s3
                    # Saving tr and ts into Directory to keep track and verify it is correct
                    Surf3['tr'].append(tr3)
                    Surf3['ts'].append(ts3)
                    # Determine if the ray impacted or not. If ts is between 0 and 1, it did impact
                    # If it impacts, then find surface with smallest tr
                    # If ts is less than 0 or bigger than 1 then no impact -> Inf
                    # Comes from Eq (9)
                    if ts3 < 0 or ts3 > 1:
                        Surf3['Impacted_surfaces'].append('Inf') # Saving non-impact rays as Inf
                    elif 0 <= ts3 <= 1:
                        Surf3['Impacted_surfaces'].append(tr3) # Saving impact rays as tr
                    # Gather Vector Data from Directory needed, so receiving rays is surface 4
                    x0_s4 = Surf4['x0,s']
                    y0_s4 = Surf4['y0,s']
                    xhat_s4 = Surf4['xhat,s']
                    yhat_s4 = Surf4['yhat,s']
                    # Calculating tr from Eq (7) by setting vector eq of surf4 and eq of rays from surf1 equals
                    tr4 = (xhat_s4 * y0_s4 - y0_ray * xhat_s4 + yhat_s4 * x0_ray - yhat_s4 * x0_s4) / (
                    yhat_ray * xhat_s4 - yhat_s4 * xhat_ray)
                    # Calculating ts from Eq (8) by setting vector eq of surf4 and eq of rays from surf1 equals
                    ts4 = (x0_ray + xhat_ray * tr4 - x0_s4) / xhat_s4
                    # Saving tr and ts into Directory to keep track and verify it is correct
                    Surf4['tr'].append(tr4)
                    Surf4['ts'].append(ts4)
                    # Determine if the ray impacted or not. If ts is between 0 and 1, it did impact
                    # If it impacts, then find surface with the smallest tr
                    # If ts is less than 0 or bigger than 1 then no impact -> Inf
                    # Comes from Eq (9)
                    if ts4 < 0 or ts4 > 1:
                        Surf4['Impacted_surfaces'].append('Inf') # Saving non-impact rays as Inf
                    elif 0 <= ts4 <= 1:
                        Surf4['Impacted_surfaces'].append(tr4) # Saving impact rays as tr
                    # Gather Vector Data from Directory needed, so receiving rays is surface 5
                    x0_s5 = Surf5['x0,s']
                    y0_s5 = Surf5['y0,s']
                    xhat_s5 = Surf5['xhat,s']
                    yhat_s5 = Surf5['yhat,s']
                    # Calculating tr from Eq (7) by setting vector eq of surf5 and eq of rays from surf1 equals
                    tr5 = (xhat_s5 * y0_s5 - y0_ray * xhat_s5 + yhat_s5 * x0_ray - yhat_s5 * x0_s5) / (
                    yhat_ray * xhat_s5 - yhat_s5 * xhat_ray)
                    # Calculating ts from Eq (8) by setting vector eq of surf5 and eq of rays from surf1 equals
                    ts5 = (x0_ray + xhat_ray * tr5 - x0_s5) / xhat_s5
                    # Saving tr and ts into Directory to keep track and verify it is correct
                    Surf5['tr'].append(tr5)
                    Surf5['ts'].append(ts5)
                    # Determine if the ray impacted or not. If ts is between 0 and 1, it did impact
                    # If it impacts, then find surface with the smallest tr
                    # If ts is less than 0 or bigger than 1 then no impact -> Inf
                    # Comes from Eq (9)
                    if ts5 < 0 or ts5 > 1:
                        Surf5['Impacted_surfaces'].append('Inf') # Saving non-impact rays as Inf
                    elif 0 <= ts5 <= 1:
                        Surf5['Impacted_surfaces'].append(tr5) # Saving impact rays as tr
                    # Gather Vector Data from Directory needed, so receiving rays is surface 6
                    x0_s6 = Surf6['x0,s']
                    y0_s6 = Surf6['y0,s']
                    xhat_s6 = Surf6['xhat,s']
                    yhat_s6 = Surf6['yhat,s']
                    # Calculating tr from Eq (7) by setting vector eq of surf6 and eq of rays from surf1 equals
                    tr6 = (xhat_s6 * y0_s6 - y0_ray * xhat_s6 + yhat_s6 * x0_ray - yhat_s6 * x0_s6) / (
                    yhat_ray * xhat_s6 - yhat_s6 * xhat_ray)
                    # Calculating ts from Eq (8) by setting vector eq of surf6 and eq of rays from surf1 equals
                    ts6 = (x0_ray + xhat_ray * tr6 - x0_s6) / xhat_s6
                    # Saving tr and ts into Directory to keep track and verify it is correct
                    Surf6['tr'].append(tr6)
                    Surf6['ts'].append(ts6)
                    # Determine if the ray impacted or not. If ts is between 0 and 1, it did impact
                    # If it impacts, then find surface with the smallest tr
                    # If ts is less than 0 or bigger than 1 then no impact -> Inf
                    # Comes from Eq (9)
                    if ts6 < 0 or ts6 > 1:
                        Surf6['Impacted_surfaces'].append('Inf') # Saving non-impact rays as Inf
                    elif 0 <= ts6 <= 1:
                        Surf6['Impacted_surfaces'].append(tr6) # Saving impact rays as tr
                    # Gather Vector Data from Directory needed, so receiving rays is surface 7
                    x0_s7 = Surf7['x0,s']
                    y0_s7 = Surf7['y0,s']
                    xhat_s7 = Surf7['xhat,s']
                    yhat_s7 = Surf7['yhat,s']
                    # Calculating tr from Eq (7) by setting vector eq of surf7 and eq of rays from surf1 equals
                    tr7 = (xhat_s7 * y0_s7 - y0_ray * xhat_s7 + yhat_s7 * x0_ray - yhat_s7 * x0_s7) / (
                    yhat_ray * xhat_s7 - yhat_s7 * xhat_ray)
                    # Calculating ts from Eq (8) by setting vector eq of surf7 and eq of rays from surf1 equals
                    ts7 = (x0_ray + xhat_ray * tr7 - x0_s7) / xhat_s7
                    # Saving tr and ts into Directory to keep track and verify it is correct
                    Surf7['tr'].append(tr7)
                    Surf7['ts'].append(ts7)
                    # Determine if the ray impacted or not. If ts is between 0 and 1, it did impact
                    # If it impacts, then find surface with the smallest tr
                    # If ts is less than 0 or bigger than 1 then no impact -> Inf
                    # Comes from Eq (9)
                    if ts7 < 0 or ts7 > 1:
                        Surf7['Impacted_surfaces'].append('Inf') # Saving non-impact rays as Inf
                    elif 0 <= ts7 <= 1:
                        Surf7['Impacted_surfaces'].append(tr7) # Saving impact rays as tr
                    # Gather Vector Data from Directory needed, so receiving rays is surface 8
                    x0_s8 = Surf8['x0,s']
                    y0_s8 = Surf8['y0,s']
                    xhat_s8 = Surf8['xhat,s']
                    yhat_s8 = Surf8['yhat,s']
                    # Calculating tr from Eq (7) by setting vector eq of surf8 and eq of rays from surf1 equals
                    tr8 = (xhat_s8 * y0_s8 - y0_ray * xhat_s8 + yhat_s8 * x0_ray - yhat_s8 * x0_s8) / (
                    yhat_ray * xhat_s8 - yhat_s8 * xhat_ray)
                    # Calculating ts from Eq (8) by setting vector eq of surf8 and eq of rays from surf1 equals
                    ts8 = (x0_ray + xhat_ray * tr8 - x0_s8) / xhat_s8
                    # Saving tr and ts into Directory to keep track and verify it is correct
                    Surf8['tr'].append(tr8)
                    Surf8['ts'].append(ts8)
                    # Determine if the ray impacted or not. If ts is between 0 and 1, it did impact
                    # If it impacts, then find surface with the smallest tr
                    # If ts is less than 0 or bigger than 1 then no impact -> Inf
                    # Comes from Eq (9)
                    if ts8 < 0 or ts8 > 1:
                        Surf8['Impacted_surfaces'].append('Inf') # Saving non-impact rays as Inf
                    elif 0 <= ts8 <= 1:
                        Surf8['Impacted_surfaces'].append(tr8) # Saving impact rays as tr
                    # Gather Vector Data from Directory needed, so receiving rays is surface 9
                    x0_s9 = Surf9['x0,s']
                    y0_s9 = Surf9['y0,s']
                    xhat_s9 = Surf9['xhat,s']
                    yhat_s9 = Surf9['yhat,s']
                    # Calculating tr from Eq (7) by setting vector eq of surf9 and eq of rays from surf1 equals
                    tr9 = (xhat_s9 * y0_s9 - y0_ray * xhat_s9 + yhat_s9 * x0_ray - yhat_s9 * x0_s9) / (
                    yhat_ray * xhat_s9 - yhat_s9 * xhat_ray)
                    # Calculating ts from Eq (8) by setting vector eq of surf9 and eq of rays from surf1 equals
                    ts9 = (x0_ray + xhat_ray * tr9 - x0_s9) / xhat_s9
                    # Saving tr and ts into Directory to keep track and verify it is correct
                    Surf9['tr'].append(tr9)
                    Surf9['ts'].append(ts9)
                    # Determine if the ray impacted or not. If ts is between 0 and 1, it did impact
                    # If it impacts, then find surface with the smallest tr
                    # If ts is less than 0 or bigger than 1 then no impact -> Inf
                    # Comes from Eq (9)
                    if ts9 < 0 or ts9 > 1:
                        Surf9['Impacted_surfaces'].append('Inf') # Saving non-impact rays as Inf
                    elif 0 <= ts9 <= 1:
                        Surf9['Impacted_surfaces'].append(tr9) # Saving impact rays as tr
                    # Gather Vector Data from Directory needed, so receiving rays is surface 10
                    x0_s10 = Surf10['x0,s']
                    y0_s10 = Surf10['y0,s']
                    xhat_s10 = Surf10['xhat,s']
                    yhat_s10 = Surf10['yhat,s']
                    # Calculating tr from Eq (7) by setting vector eq of surf10 and eq of rays from surf1 equals
                    tr10 = (xhat_s10 * y0_s10 - y0_ray * xhat_s10 + yhat_s10 * x0_ray - yhat_s10 * x0_s10) / (
                    yhat_ray * xhat_s10 - yhat_s10 * xhat_ray)
                    # Calculating ts from Eq (8) by setting vector eq of surf10 and eq of rays from surf1 equals
                    ts10 = (x0_ray + xhat_ray * tr10 - x0_s10) / xhat_s10
                    # Saving tr and ts into Directory to keep track and verify it is correct
                    Surf10['tr'].append(tr10)
                    Surf10['ts'].append(ts10)
                    # Determine if the ray impacted or not. If ts is between 0 and 1, it did impact
                    # If it impacts, then find surface with the smallest tr
                    # If ts is less than 0 or bigger than 1 then no impact -> Inf
                    # Comes from Eq (9)
                    if ts10 < 0 or ts10 > 1:
                        Surf10['Impacted_surfaces'].append('Inf') # Saving non-impact rays as Inf
                    elif 0 <= ts10 <= 1:
                        Surf10['Impacted_surfaces'].append(tr10) # Saving impact rays as tr
                    # Gather Vector Data from Directory needed, so receiving rays is surface 11
                    x0_s11 = Surf11['x0,s']
                    y0_s11 = Surf11['y0,s']
                    xhat_s11 = Surf11['xhat,s']
                    yhat_s11 = Surf11['yhat,s']
                    # Calculating tr from Eq (7) by setting vector eq of surf11 and eq of rays from surf1 equals
                    tr11 = (xhat_s11 * y0_s11 - y0_ray * xhat_s11 + yhat_s11 * x0_ray - yhat_s11 * x0_s11) / (
                    yhat_ray * xhat_s11 - yhat_s11 * xhat_ray)
                    # Calculating ts from Eq (8) by setting vector eq of surf11 and eq of rays from surf1 equals
                    ts11 = (x0_ray + xhat_ray * tr11 - x0_s11) / xhat_s11
                    # Saving tr and ts into Directory to keep track and verify it is correct
                    Surf11['tr'].append(tr11)
                    Surf11['ts'].append(ts11)
                    # Determine if the ray impacted or not. If ts is between 0 and 1, it did impact
                    # If it impacts, then find surface with the smallest tr
                    # If ts is less than 0 or bigger than 1 then no impact -> Inf
                    # Comes from Eq (9)
                    if ts11 < 0 or ts11 > 1:
                        Surf11['Impacted_surfaces'].append('Inf') # Saving non-impact rays as Inf
                    elif 0 <= ts11 <= 1:
                        Surf11['Impacted_surfaces'].append(tr11) # Saving impact rays as tr
                        
        # Calculating Total number of rayts emitted (should be N_rays as well)
        total_rays = len(Surf2['Impacted_surfaces'])
        
        # # Calculating Number of Impacted rays by counting entries in Impacted_Surfaces that are not Inf
        # # meaning the amount of rays that hit a surface
        # num_impacted_rays_2 = sum(1 for surface in Surf2['Impacted_surfaces'] if surface != 'Inf')
        # num_impacted_rays_3 = sum(1 for surface in Surf3['Impacted_surfaces'] if surface != 'Inf')
        # num_impacted_rays_4 = sum(1 for surface in Surf4['Impacted_surfaces'] if surface != 'Inf')
        # num_impacted_rays_5 = sum(1 for surface in Surf5['Impacted_surfaces'] if surface != 'Inf')
        # num_impacted_rays_6 = sum(1 for surface in Surf6['Impacted_surfaces'] if surface != 'Inf')
        # num_impacted_rays_7 = sum(1 for surface in Surf7['Impacted_surfaces'] if surface != 'Inf')
        # num_impacted_rays_8 = sum(1 for surface in Surf8['Impacted_surfaces'] if surface != 'Inf')
        # num_impacted_rays_9 = sum(1 for surface in Surf9['Impacted_surfaces'] if surface != 'Inf')
        # num_impacted_rays_10 = sum(1 for surface in Surf10['Impacted_surfaces'] if surface != 'Inf')
        # num_impacted_rays_11 = sum(1 for surface in Surf11['Impacted_surfaces'] if surface != 'Inf')
        
        # # Calculate the view factor F from Eq(10) and save into list to plot
        # view_factor_1_2 = num_impacted_rays_2 / total_rays
        # view_factor_1_3 = num_impacted_rays_3 / total_rays
        # view_factor_1_4 = num_impacted_rays_4 / total_rays
        # view_factor_1_5 = num_impacted_rays_5 / total_rays
        # view_factor_1_6 = num_impacted_rays_6 / total_rays
        # view_factor_1_7 = num_impacted_rays_7 / total_rays
        # view_factor_1_8 = num_impacted_rays_8 / total_rays
        # view_factor_1_9 = num_impacted_rays_9 / total_rays
        # view_factor_1_10 = num_impacted_rays_10 / total_rays
        # view_factor_1_11 = num_impacted_rays_11 / total_rays
        # view_factor_values_1_2.append(view_factor_1_2)
        # view_factor_values_1_3.append(view_factor_1_3)
        # view_factor_values_1_4.append(view_factor_1_4)
        # view_factor_values_1_5.append(view_factor_1_5)
        # view_factor_values_1_6.append(view_factor_1_6)
        # view_factor_values_1_7.append(view_factor_1_7)
        # view_factor_values_1_8.append(view_factor_1_8)
        # view_factor_values_1_9.append(view_factor_1_9)
        # view_factor_values_1_10.append(view_factor_1_10)
        # view_factor_values_1_11.append(view_factor_1_11)
        
        # Calculating Number of Impacted rays by counting entries in Impacted_Surfaces that are not Inf
        num_impacted_rays = {i: sum(1 for surface in surf['Impacted_surfaces'] if surface != 'Inf') for i, surf in enumerate(
           [Surf2, Surf3, Surf4, Surf5, Surf6, Surf7, Surf8, Surf9, Surf10, Surf11], start=2)}
       
       # Calculate the view factor F from Eq(10) and save into list to plot
        view_factors = {i: num_impacted / total_rays for i, num_impacted in num_impacted_rays.items()}
        
        # Add the results to rows
        row = {'V_groove_angle': np.rad2deg(V_groove_angle_half * 2)}
        row.update(view_factors)
        rows.append(row)

# Create DataFrame and save to Excel
df = pd.DataFrame(rows)
df.to_excel('view_factors.xlsx', index=False)