import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import radians, degrees

# Constants
LATITUDE = 39.7589
LONGITUDE = -84.1916
GT = 0.95
A = 1
ETA = 0.95
H = 10
TC_OFFSET = 1
ALBEDO = 0.2  # Albedo factor

# Tracking parameters
TRACKING_ANGLE = 0.0  # Initial tracking angle (degrees)
TRACKING_FACTOR = 0.5  # Tracking factor (0.0 for fixed, 1.0 for full tracking)

# Read data from CSV file
def read_data(file_name):
    data = pd.read_csv(file_name)
    return data

# Calculate solar angles with tracking
def calculate_solar_angles(data):
    # Extract year, month, day, and hour columns and combine them into a single datetime column
    data['Datetime'] = pd.to_datetime(data[['Year', 'Month', 'Day', 'Hour']])
    
    # Calculate solar angles for the entire year
    data['DayOfYear'] = data['Datetime'].dt.dayofyear
    data['HourOfYear'] = (data['DayOfYear'] - 1) * 24 + data['Hour']
    
    # Calculate solar declination
    data['Declination'] = np.radians(data['DayOfYear'].apply(lambda doy: 23.45 * np.sin(np.radians(360 / 365 * (doy - 81)))))
    
    # Calculate solar hour angle with tracking
    data['SolarHourAngle'] = 15 * (data['HourOfYear'] - 12)
    
    # Calculate solar elevation angle
    sin_elevation_angle = np.sin(np.radians(LATITUDE)) * np.sin(data['Declination']) + np.cos(np.radians(LATITUDE)) * np.cos(data['Declination']) * np.cos(np.radians(data['SolarHourAngle']))
    data['SolarElevationAngle'] = np.degrees(np.arcsin(sin_elevation_angle))
    
    # Calculate solar azimuth angle
    sin_azimuth_angle = -np.cos(np.radians(data['SolarHourAngle']))
    cos_azimuth_angle = -np.sin(np.radians(LATITUDE)) * np.cos(data['Declination']) + np.cos(np.radians(LATITUDE)) * np.sin(data['Declination']) * np.cos(np.radians(data['SolarHourAngle']))
    data['SolarAzimuthAngle'] = np.degrees(np.arctan2(sin_azimuth_angle, cos_azimuth_angle))
    
    return data

# Calculate solar irradiance
def calculate_solar_irradiance(data):
    data['Irradiance'] = GT * data['GHI (W/m^2)']
    return data

# Calculate power output
def calculate_power_output(data):
    data['PowerOutput'] = A * data['Irradiance'] * ETA * np.cos(np.radians(data['SolarElevationAngle'])) * (1 - TRACKING_FACTOR * np.abs(np.cos(np.radians(data['SolarAzimuthAngle'] - TRACKING_ANGLE)))) * (1 - ALBEDO)
    return data

# Plot power output by month
def plot_power_output_by_month(data):
    # Group the data by month and plot each month separately
    grouped = data.groupby(data['Datetime'].dt.month)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    for month, group in grouped:
        ax.plot(group['Datetime'], group['PowerOutput'], label=f'Month {month}')
    
    plt.xlabel('Time')
    plt.ylabel('Power Output (W)')
    plt.title('Hourly DC Power Production by Month with Single-Axis Tracking and Albedo Factor')
    plt.legend(title='Months', labels=['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'])
    plt.grid(True)
    plt.show()

# Save data to CSV file
def save_data(data, file_name):
    save_path = 'D:\Graduate\gpt-engineer\solar_app\workspace\\' + file_name
    data.to_csv(save_path)

# Main function
def main():
    data = read_data('Dayton_tmy_2021.csv')
    data = calculate_solar_angles(data)
    data = calculate_solar_irradiance(data)
    data = calculate_power_output(data)
    plot_power_output_by_month(data)
    save_data(data, 'SolarAPP_Tracking_Albedo.csv')

if __name__ == "__main__":
    main()
