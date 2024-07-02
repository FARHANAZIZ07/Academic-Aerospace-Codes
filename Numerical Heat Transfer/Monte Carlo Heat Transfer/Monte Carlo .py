# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 14:55:30 2023

@author: Farhan
"""

import random
import matplotlib.pyplot as plt
import numpy as np

radius = 1

def generate_random_point():
    x = random.uniform (0, radius)
    y = random.uniform(0, radius)
    return x, y

for max_dots in [1000, 10000]:
    count_points_inside = 0
    count_points_outside = 0
    points_inside = []
    points_outside= []
    pi_values = []
   
    for i in range(max_dots):
        x, y = generate_random_point()
        if x**2 + y**2 <= radius**2: # x^2 + y^2 = r^2
            points_inside.append((x, y))
            count_points_inside += 1
        else:
            points_outside.append((x, y))
            count_points_outside += 1
            
        current_total_points = count_points_inside + count_points_outside
        pi = 4 * count_points_inside / current_total_points
        pi_values.append(pi)
    
    print (f"Final Estimate: pi = {pi} for {max_dots} points")

# To generate radius line
    theta = np. linspace (0, 2*np.pi, 100)
    x_radius = radius * np.cos(theta)
    y_radius = radius * np.sin(theta)
    points = np.linspace(0, max_dots+1,max_dots)
    
    #points inside and outside of the radius plot
    if max_dots == 1000:
        plt.figure(figsize=(6, 6))
        plt.scatter (*zip(*points_inside), color='blue', linewidth=0.2, label='Inside of the Circle') 
        plt.scatter (*zip(*points_outside), color='red', linewidth=0.2, label='Outside of the Circle') 
        plt.plot(x_radius, y_radius, linewidth=1.7, color='black', label='Circle')
        plt. legend (loc= 'upper left', bbox_to_anchor=(1, 1))
        plt.title(f'Monte Carlo Method to Calculate Pi for {max_dots} dots')
        plt.xlabel('X-axis (radius)')
        plt.ylabel('Y-axis (radius)') 
        plt.xlim (0, radius) 
        plt.ylim(0, radius)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()
        #pi values throughout the max dots
    if max_dots == 10000:
        plt.figure()
        plt.plot(points,pi_values, linewidth=1.7, color='black', label='Radius')
        plt.legend()
        plt.title(f'Pi values as function of number of points for {max_dots} dots') 
        plt.xlabel('Number of points')
        plt.ylabel('Pi value')
        plt.xlim(0, max_dots+1)
        plt.ylim (0, 4)
        plt.show()