# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 21:27:03 2024

@author: Farhan
"""

import pandas as pd
import numpy as np
import copy
import matplotlib.pyplot as plt
# Load data from the Excel file with view factors
excel_file_path = r"C:\Users\............\view_factors.xlsx"

# Load data from different sheets
Numerical_view_factors = pd.read_excel(excel_file_path, sheet_name="Numerical View Factors")

# Extract Data for view factors for different V groove angles
View_factors_30_degrees = Numerical_view_factors['30 degrees'].tolist()
View_factors_90_degrees = Numerical_view_factors['90 degrees'].tolist()
View_factors_120_degrees = Numerical_view_factors['120 degrees'].tolist()
View_factors_30_degrees_surface1 = View_factors_30_degrees[0:10]
View_factors_90_degrees_surface1 = View_factors_90_degrees[0:10]
View_factors_120_degrees_surface1 = View_factors_120_degrees[0:10]
View_factors_30_degrees_surface2 = View_factors_30_degrees[10:20]
View_factors_90_degrees_surface2 = View_factors_90_degrees[10:20]
View_factors_120_degrees_surface2 = View_factors_120_degrees[10:20]
View_factors_30_degrees_surface3 = View_factors_30_degrees[20:30]
View_factors_90_degrees_surface3 = View_factors_90_degrees[20:30]
View_factors_120_degrees_surface3 = View_factors_120_degrees[20:30]
View_factors_30_degrees_surface4 = View_factors_30_degrees[30:40]
View_factors_90_degrees_surface4 = View_factors_90_degrees[30:40]
View_factors_120_degrees_surface4 = View_factors_120_degrees[30:40]
View_factors_30_degrees_surface5 = View_factors_30_degrees[40:50]
View_factors_90_degrees_surface5 = View_factors_90_degrees[40:50]
View_factors_120_degrees_surface5 = View_factors_120_degrees[40:50]
View_factors_30_degrees_surface6 = View_factors_30_degrees[50:60]
View_factors_90_degrees_surface6 = View_factors_90_degrees[50:60]
View_factors_120_degrees_surface6 = View_factors_120_degrees[50:60]
View_factors_30_degrees_surface7 = View_factors_30_degrees[60:70]
View_factors_90_degrees_surface7 = View_factors_90_degrees[60:70]
View_factors_120_degrees_surface7 = View_factors_120_degrees[60:70]
View_factors_30_degrees_surface8 = View_factors_30_degrees[70:80]
View_factors_90_degrees_surface8 = View_factors_90_degrees[70:80]
View_factors_120_degrees_surface8 = View_factors_120_degrees[70:80]
View_factors_30_degrees_surface9 = View_factors_30_degrees[80:90]
View_factors_90_degrees_surface9 = View_factors_90_degrees[80:90]
View_factors_120_degrees_surface9 = View_factors_120_degrees[80:90]
View_factors_30_degrees_surface10 = View_factors_30_degrees[90:100]
View_factors_90_degrees_surface10 = View_factors_90_degrees[90:100]
View_factors_120_degrees_surface10 = View_factors_120_degrees[90:100]
View_factors_30_degrees_by_surfaces = [View_factors_30_degrees_surface1, View_factors_30_degrees_surface2,
View_factors_30_degrees_surface3 ,View_factors_30_degrees_surface4, View_factors_30_degrees_surface5, View_factors_30_degrees_surface6,
View_factors_30_degrees_surface7, View_factors_30_degrees_surface8, View_factors_30_degrees_surface9, View_factors_30_degrees_surface10]
View_factors_90_degrees_by_surfaces = [View_factors_90_degrees_surface1, View_factors_90_degrees_surface2,
View_factors_90_degrees_surface3 ,View_factors_90_degrees_surface4, View_factors_90_degrees_surface5, View_factors_90_degrees_surface6,
View_factors_90_degrees_surface7, View_factors_90_degrees_surface8, View_factors_90_degrees_surface9, View_factors_90_degrees_surface10]
View_factors_120_degrees_by_surfaces = [View_factors_120_degrees_surface1, View_factors_120_degrees_surface2,
View_factors_120_degrees_surface3 ,View_factors_120_degrees_surface4, View_factors_120_degrees_surface5,
View_factors_120_degrees_surface6, View_factors_120_degrees_surface7, View_factors_120_degrees_surface8,
View_factors_120_degrees_surface9, View_factors_120_degrees_surface10]
View_factors_degrees = [View_factors_30_degrees_by_surfaces, View_factors_90_degrees_by_surfaces, View_factors_120_degrees_by_surfaces]
#knowns from Radiosity balance equation
epsilon_values = [0.1, 0.5, 0.9] #emissivity values
sigma = 5.67e-8 #stefan boltzmann constant
T_i = 300 #Temperature of groove walls in Kelvin
L_surface = 0.1 # length of groove wall in m
N_surf = 10 #number of surfaces in a wall
w = L_surface/N_surf # lenght of each surface segment in m
A_i = w #area of each surface, in this case since it is only w for 1D
#For Gauss Seidel
Convergence_criteria = 0.000000000001
iteration = 0
J_final_all = [] #To store final Radiosity for different angles and emissivity values
for View_factors_x_degrees_by_surfaces in View_factors_degrees: #For different angles
    for epsilon in epsilon_values: #For different emissivity values
        J_guess = [[1 for j in range (0,N_surf)] for i in range(0,N_surf)] #initial guess
        converged = False
        J_final = [] #To store values
        #Gauss-Seidel Iteration
        while converged == False:
            J_new_values = [] #To store values
            J_new_list = [] #To store values
            J_new = copy.deepcopy(J_guess) #Set J_new to be last guess
            sum_previous_term_values = [] #To store values
            for i in range(N_surf):
                previous_term = [] #To store values
                for j in range(N_surf):
                    term_ij = J_guess[i][j]*View_factors_x_degrees_by_surfaces[i][j]*A_i #Find term inside summation from previous guess
                    previous_term.append(term_ij) #To store values
                    sum_previous_term = sum(previous_term) #Summation of all surfaces
                    sum_previous_term_values.append(sum_previous_term) #To store values
                    
                J_new = (epsilon*A_i*sigma*T_i**4+(1-epsilon)*sum_previous_term)/A_i #Calculating the new radiosity
                J_new_values.append(J_new) #To store values
                J_new_list.append([J_new]*10) #To store values
                
            iteration = iteration+1
            converged = True
            J_new_array = np.array(J_new_list)
            J_guess_array = np.array(J_guess)
            
            #Check if converged or not
            for i in range(0,N_surf):
                for j in range(0,N_surf):
                    if abs(abs(J_new_array[i][j])-abs(J_guess_array[i][j])) > Convergence_criteria:
                        converged = False
            J_guess = copy.deepcopy(J_new_list) #Saving new J_guess as the last J_new
        for entry in J_new_list: #Save final radiosity as lists for each angle and emissivity values
            J_final.append(entry[0])
        J_final_all.append(J_final)
        #Index 0-2 for angle = 30 (Index 0 for epsilon 0.1, Index 1 for epsilon 0.5, Index 2 for epsilon 0.9)
        #Index 3-5 for angle = 90 (Index 3 for epsilon 0.1, Index 4 for epsilon 0.5, Index 5 for epsilon 0.9)
        #Index 6-8 for angle = 120 (Index 6 for epsilon 0.1, Index 7 for epsilon 0.5, Index 8 for epsilon 0.9)

#Plotting in the middle of each surface
wall_x = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5]

# Plot for 30 degrees
plt.figure(figsize=(10, 6))
plt.plot(wall_x, J_final_all[0], label='Emissivity value = 0.1')
plt.plot(wall_x, J_final_all[1], label='Emissivity value = 0.5')
plt.plot(wall_x, J_final_all[2], label='Emissivity value = 0.9')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xlabel('Wall Position (m)')
plt.ylabel('Radiosity J (W/m2)')
plt.title('Radiosity vs Wall Position for a V-groove angle of 30 degrees')
# Plot for 90 degrees
plt.figure(figsize=(10, 6))
plt.plot(wall_x, J_final_all[3], label='Emissivity value = 0.1')
plt.plot(wall_x, J_final_all[4], label='Emissivity value = 0.5')
plt.plot(wall_x, J_final_all[5], label='Emissivity value = 0.9')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xlabel('Wall Position (m)')
plt.ylabel('Radiosity J (W/m2)')
plt.title('Radiosity vs Wall Position for a V-groove angle of 90 degrees')
# Plot for 120 degrees
plt.figure(figsize=(10, 6))
plt.plot(wall_x, J_final_all[6], label='Emissivity value = 0.1')
plt.plot(wall_x, J_final_all[7], label='Emissivity value = 0.5')
plt.plot(wall_x, J_final_all[8], label='Emissivity value = 0.9')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xlabel('Wall Position (m)')
plt.ylabel('Radiosity J (W/m2)')
plt.title('Radiosity vs Wall Position for a V-groove angle of 120 degrees')
#Plot all together
plt.figure(figsize=(10, 6))
plt.plot(wall_x, J_final_all[0], label='Angle = 30 Degrees, Emissivity value = 0.1')
plt.plot(wall_x, J_final_all[1], label='Angle = 30 Degrees, Emissivity value = 0.5')
plt.plot(wall_x, J_final_all[2], label='Angle = 30 Degrees, Emissivity value = 0.9')
plt.plot(wall_x, J_final_all[3], label='Angle = 90 Degrees, Emissivity value = 0.1')
plt.plot(wall_x, J_final_all[4], label='Angle = 90 Degrees, Emissivity value = 0.5')
plt.plot(wall_x, J_final_all[5], label='Angle = 90 Degrees, Emissivity value = 0.9')
plt.plot(wall_x, J_final_all[6], label='Angle = 120 Degrees, Emissivity value = 0.1')
plt.plot(wall_x, J_final_all[7], label='Angle = 120 Degrees, Emissivity value = 0.5')
plt.plot(wall_x, J_final_all[8], label='Angle = 120 Degrees, Emissivity value = 0.9')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xlabel('Wall Position (m)')
plt.ylabel('Radiosity J (W/m2)')
plt.title('Radiosity vs Wall Position for various V-groove angles and emissivities')
plt.show()

# For extra credit:
J_final_30 = [J_final_all[0], J_final_all[1], J_final_all[2]]
J_final_90 = [J_final_all[3], J_final_all[4], J_final_all[5]]
J_final_120 = [J_final_all[6], J_final_all[7], J_final_all[8]]
J_i = [J_final_30, J_final_90, J_final_120] #Final emissivities, first is for angles 30, 90, 120, then inside for emissivities 0.1, 0.5,
0.9
View_factors_leaving_top_surface_values = [] #To store values

#Sum all view factors on each surface and subtract by 1 to find the view factor from each surface on left side to top surface from summation
for View_factors_x_degrees_by_surfaces in View_factors_degrees: #For different angles
    for i in View_factors_x_degrees_by_surfaces:
        View_factor_leaving_top_surface = 1-sum(i)
        View_factors_leaving_top_surface_values.append(View_factor_leaving_top_surface)

#Saving all view factors for easier use
V_F_30 = [View_factors_leaving_top_surface_values[0], View_factors_leaving_top_surface_values[1],
View_factors_leaving_top_surface_values[2],View_factors_leaving_top_surface_values[3],View_factors_leaving_top_surface_values[4],
View_factors_leaving_top_surface_values[5],View_factors_leaving_top_surface_values[6],View_factors_leaving_top_surface_values[7],
View_factors_leaving_top_surface_values[8],View_factors_leaving_top_surface_values[9]]
V_F_90 = [View_factors_leaving_top_surface_values[10], View_factors_leaving_top_surface_values[11],
View_factors_leaving_top_surface_values[12],View_factors_leaving_top_surface_values[13],View_factors_leaving_top_surface_values[14],
View_factors_leaving_top_surface_values[15],View_factors_leaving_top_surface_values[16],View_factors_leaving_top_surface_values[17],
View_factors_leaving_top_surface_values[18],View_factors_leaving_top_surface_values[19]]
V_F_120 = [View_factors_leaving_top_surface_values[20], View_factors_leaving_top_surface_values[21],
View_factors_leaving_top_surface_values[22],View_factors_leaving_top_surface_values[23],View_factors_leaving_top_surface_values[24],
View_factors_leaving_top_surface_values[25],View_factors_leaving_top_surface_values[26],View_factors_leaving_top_surface_values[27],
View_factors_leaving_top_surface_values[28],View_factors_leaving_top_surface_values[29]]
View_factors_leaving_list_30 = [V_F_30 ,V_F_30, V_F_30]
View_factors_leaving_list_90 = [V_F_90 ,V_F_90, V_F_90]
View_factors_leaving_list_120 = [V_F_120 ,V_F_120, V_F_120]
View_factors_leaving_lists = [View_factors_leaving_list_30 ,View_factors_leaving_list_90, View_factors_leaving_list_120]
F_i_top = View_factors_leaving_lists #View factors from surface x to top surface for each angle first then for each emissivity (all same
# bc view factor only affected from geometry)
    
#Calculating total heat loss for different angles and emissivities
q_loss_for_every_angle = []
q_total_loss = []

for i in range(len(View_factors_degrees)):
    q_loss_for_every_epsilon= []
    q_total_loss_per_surfaces = []
    for j in range(len(epsilon_values)):
        q_loss_per_surfaces = []
        for k in range(0,10):
            q_loss_per_surface = J_i[i][j][k]*View_factors_leaving_lists[i][j][k]*w
            q_loss_per_surfaces.append(q_loss_per_surface) #Heat loss in every surface
        q_loss_for_every_epsilon.append(q_loss_per_surfaces)
        q_total_loss_per_surface = sum(q_loss_per_surfaces) #Adding all the heat loss from every surface
        q_total_loss_per_surfaces.append(q_total_loss_per_surface)
    q_loss_for_every_angle.append(q_loss_for_every_epsilon)
    q_total_loss.append(q_total_loss_per_surfaces) #total heat loss first for angles 30, 90, 120 and then inside for epsilon 0.1, 0.5, 0.9

q_total_loss_epsilon_1 = []
q_total_loss_epsilon_2 = []
q_total_loss_epsilon_3 = []
q_total_loss_epsilon_1_good = []
q_total_loss_epsilon_2_good = []
q_total_loss_epsilon_3_good = []
for entry in q_total_loss:
    q_total_loss_epsilon_1.append(entry[0])
    q_total_loss_epsilon_2.append(entry[1])
    q_total_loss_epsilon_3.append(entry[2])
    
#Plotting total heat loss for different angles and emissivities
angles = [30, 90, 120]
# Plot
plt.figure(figsize=(10, 6))
plt.scatter(angles, q_total_loss_epsilon_1, label='Emissivity value = 0.1', marker='o')
plt.scatter(angles, q_total_loss_epsilon_2, label='Emissivity value = 0.5', marker='o')
plt.scatter(angles, q_total_loss_epsilon_3, label='Emissivity value = 0.9', marker='o')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xlabel('Angles (Degrees)')
plt.ylabel('Total Heat loss (W/m) ')
plt.title('Total Heat loss vs V-groove Angle for different emissivities')
plt.show()