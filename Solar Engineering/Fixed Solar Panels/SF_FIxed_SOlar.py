# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 03:07:44 2023

@author: Farhan
"""

import pandas as pd
import pvlib

# Define the layout

# Read input values
weather_data = pd.read_csv("Dayton_tmy_2021.csv")
panel_efficiency = 0.19
inverter_efficiency = 0.96
area = 1.6
albedo = 0.2
surface_tilt = 30
surface_azimuth = 180
latitude = 39.417
longitude = -84.17

# Read weather data

time_data = weather_data[['Year', 'Month', 'Day', 'Hour', 'Minute']].astype(str)
time_range = pd.to_datetime(time_data['Year'] + '-' + time_data['Month'] + '-' + time_data['Day'] + ' ' + time_data['Hour'] + ':' + time_data['Minute'])
time_range = time_range.reset_index(drop=True)  # Reset index for alignment

ghi = weather_data['GHI'].reset_index(drop=True)  # Reset index for alignment
dni = weather_data['DNI'].reset_index(drop=True)  # Reset index for alignment
dhi = weather_data['DHI'].reset_index(drop=True)  # Reset index for alignment

# Solar position
solar_position = pvlib.solarposition.get_solarposition(time_range, latitude, longitude)
solar_position = solar_position.reset_index(drop=True)  # Reset index for alignment

aoi = pvlib.irradiance.aoi(surface_tilt, surface_azimuth, solar_position['apparent_zenith'], solar_position['azimuth'])

# Add AOI to the solar_position DataFrame
solar_position['aoi'] = aoi

# Calculate POA irradiance
poa_irradiance = pvlib.irradiance.get_total_irradiance(surface_tilt, surface_azimuth, solar_position['apparent_zenith'], solar_position['azimuth'], dni, ghi, dhi, albedo=albedo)['poa_global']

# Calculate DC output
dc_output = poa_irradiance * area * panel_efficiency

# Calculate AC output
ac_output = dc_output * inverter_efficiency

# Save DC and AC output to CSV
output_df = pd.DataFrame({
    'DC_Output': dc_output,
    'AC_Output': ac_output
})
output_df.to_csv('dc_ac_output.csv', index=False)

# Calculate annual energy production and peak power for DC
annual_energy_production_dc = dc_output.sum() / 1000  # Convert to kWh
peak_power_dc = dc_output.max()  # in kW

# Calculate annual energy production and peak power for AC
annual_energy_production_ac = ac_output.sum() / 1000  # Convert to kWh
peak_power_ac = ac_output.max()  # in kW

# Print to console
print(f'Annual Energy Production (DC): {annual_energy_production_dc} kWh')
print(f'Peak Power (DC): {peak_power_dc} kW')
print(f'Annual Energy Production (AC): {annual_energy_production_ac} kWh')
print(f'Peak Power (AC): {peak_power_ac} kW')

# Show in popup
print(f'Annual Energy Production (DC): {annual_energy_production_dc} kWh', f'Peak Power (DC): {peak_power_dc} kW',
         f'Annual Energy Production (AC): {annual_energy_production_ac} kWh', f'Peak Power (AC): {peak_power_ac} kW')
