import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pvlib

# Define parameters for the shading scenario
num_rows = 2  # Number of rows of panels
num_panels_per_row = 3  # Number of panels per row
roof_height = 40  # Height of the roof (in feet)
tree_height = 60  # Height of the tree (in feet)
tree_distance = 20  # Distance of the tree from the house (in feet)
panel_length = 6  # Length of each panel (in feet)
panel_width = 3  # Width of each panel (in feet)
roof_pitch_angle = 30  # Roof pitch angle in degrees

# Create time intervals (e.g., hourly data for a year)
time_intervals = pd.date_range(start='2023-01-01', end='2023-12-31', freq='H')  # Hourly data for a year

# Dayton, Ohio latitude and longitude
latitude = 39.7589  # Latitude for Dayton, Ohio
longitude = -84.1916  # Longitude for Dayton, Ohio

# Create PV system and location for PVLib calculations
system = pvlib.pvsystem.PVSystem(surface_tilt=roof_pitch_angle, surface_azimuth=180)
location = pvlib.location.Location(latitude=latitude, longitude=longitude)

# Calculate solar position for the entire year
solar_position = location.get_solarposition(time_intervals)

# Initialize arrays to store shading percentages over time
shading_percentage_over_time_upper = np.zeros(len(time_intervals))
shading_percentage_over_time_lower = np.zeros(len(time_intervals))

# Define the leaf state: 0 - no leaves, 1 - partial leaves, 2 - full leaves
leaf_state = np.zeros(len(time_intervals))

# Define seasonal leaf changes
leaf_state[1200:1700] = 2  # Summer: May to August
leaf_state[1700:2200] = 1  # Autumn: September to November
leaf_state[2200:2900] = 0  # Winter: December to March
leaf_state[2900:3000] = 1  # Spring: April

# Calculate shading percentage for upper and lower panels over time
for i in range(len(time_intervals)):
    # Calculate distances to the tree for upper and lower panels at this time interval
    if i < 1700 or (i >= 2200 and i < 2900):
        tree_distance_upper = 20  # Distance when tree has full leaves or partial leaves
    else:
        tree_distance_upper = 25  # Distance when tree has no leaves

    distance_to_tree_upper = np.sqrt(tree_distance_upper ** 2 + (roof_height - tree_height) ** 2)
    distance_to_tree_lower = np.sqrt(tree_distance ** 2 + (roof_height + 3 - tree_height) ** 2)

    # Calculate shading percentage for upper and lower panels at this time interval
    shading_percentage_upper = max(0, min(1, (panel_length / 2) / distance_to_tree_upper))
    shading_percentage_lower = max(0, min(1, (panel_length / 2) / distance_to_tree_lower))

    # Adjust shading percentage based on leaf state
    if leaf_state[i] == 1:
        shading_percentage_upper *= 0.7  # Partial leaves reduce shading
        shading_percentage_lower *= 0.7
    elif leaf_state[i] == 0:
        shading_percentage_upper *= 0.4  # No leaves significantly reduce shading
        shading_percentage_lower *= 0.4

    # Store shading percentages in the arrays
    shading_percentage_over_time_upper[i] = shading_percentage_upper
    shading_percentage_over_time_lower[i] = shading_percentage_lower

# Define inverter efficiencies
inverter_efficiency_traditional = 0.96
inverter_efficiency_string = 0.98
inverter_efficiency_micro = 0.99

# Calculate DC power generation for each inverter type
dc_output_traditional = (1 - shading_percentage_over_time_upper) * (1 - shading_percentage_over_time_lower)
dc_output_string = (1 - shading_percentage_over_time_upper) * (1 - shading_percentage_over_time_lower)
dc_output_micro = (1 - shading_percentage_over_time_upper) * (1 - shading_percentage_over_time_lower)

# Calculate AC power generation for each inverter type
ac_output_traditional = dc_output_traditional * inverter_efficiency_traditional
ac_output_string = dc_output_string * inverter_efficiency_string
ac_output_micro = dc_output_micro * inverter_efficiency_micro

# Create a DataFrame to store the data
data = pd.DataFrame({
    'Time (hours)': time_intervals,
    'Shading Percentage (Upper Panels)': shading_percentage_over_time_upper,
    'Shading Percentage (Lower Panels)': shading_percentage_over_time_lower,
    'AC Power Generation (Traditional Inverter) (kW)': ac_output_traditional,
    'AC Power Generation (String Inverter) (kW)': ac_output_string,
    'AC Power Generation (Micro-Inverter) (kW)': ac_output_micro
})

# Calculate monthly data
monthly_data = data.set_index('Time (hours)').resample('M').sum()

# Print annual energy generation for each inverter type
annual_energy_traditional = np.sum(ac_output_traditional)  # kWh
annual_energy_string = np.sum(ac_output_string)  # kWh
annual_energy_micro = np.sum(ac_output_micro)  # kWh

print(f'Annual Energy Generation (Traditional Inverter): {annual_energy_traditional} kWh')
print(f'Annual Energy Generation (String Inverter): {annual_energy_string} kWh')
print(f'Annual Energy Generation (Micro-Inverter): {annual_energy_micro} kWh')

# Save the data to a CSV file (append mode)
data.to_csv('Solar_single_axis.csv', mode='a', header=False, index=False)

# Generate plots for shading percentage over time
plt.figure(figsize=(12, 6))
plt.plot(time_intervals, shading_percentage_over_time_upper, label='Upper Panels')
plt.plot(time_intervals, shading_percentage_over_time_lower, label='Lower Panels')
plt.xlabel('Time (hours)')
plt.ylabel('Shading Percentage')
plt.title('Shading Percentage Over Time')
plt.legend()
plt.grid(True)
plt.show()

# Generate plots for total power generation for each inverter type
plt.figure(figsize=(12, 6))
plt.plot(time_intervals, ac_output_traditional, label='Traditional Inverter')
plt.plot(time_intervals, ac_output_string, label='String Inverter')
plt.plot(time_intervals, ac_output_micro, label='Micro-Inverter')
plt.xlabel('Time (hours)')
plt.ylabel('AC Power Generation (kW)')
plt.title('Total Power Generation for Different Inverter Types')
plt.legend()
plt.grid(True)
plt.show()

# Generate plots for monthly total power generation for each inverter type
plt.figure(figsize=(12, 6))
plt.plot(monthly_data.index, monthly_data['AC Power Generation (Traditional Inverter) (kW)'], label='Traditional Inverter')
plt.plot(monthly_data.index, monthly_data['AC Power Generation (String Inverter) (kW)'], label='String Inverter')
plt.plot(monthly_data.index, monthly_data['AC Power Generation (Micro-Inverter) (kW)'], label='Micro-Inverter')
plt.xlabel('Month')
plt.ylabel('Monthly Total AC Power Generation (kWh)')
plt.title('Monthly Total Power Generation for Different Inverter Types')
plt.xticks(rotation=45)
plt.legend()
plt.grid(True)
plt.show()
