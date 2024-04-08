# This is a script for creating the simulated dugong sightings test data from
# the original dugong sightings data. The goal is to create a test example dataset that
# we can make available with this repository without releasing the original
# data.
import pandas as pd
import numpy as np
from datetime import timedelta

# Load the original dataset
original_data_path = '../../original/TIDES_input-data-for-eReef_nGBR2023.csv'
data = pd.read_csv(original_data_path)

# Select a random subset of about 10% of the rows
subset_percentage = 0.1  # For example, select 10% of the data
data = data.sample(frac=subset_percentage).reset_index(drop=True)

# Apply a single random shift to all dates
date_shift = np.random.randint(-30, 31)  # Random days shift
data['date'] = pd.to_datetime(data['date'], format='%d/%m/%Y') + timedelta(days=date_shift)

# Apply individual random shifts to latitude and longitude
data['lat'] += np.random.uniform(-0.01, 0.01, size=len(data))
data['long'] += np.random.uniform(-0.02, 0.02, size=len(data))

def offset_time(group):
    # Calculate the offset once for each group
    offset = np.random.randint(-30, 31)  # Random minutes shift
    # Apply the offset to each row in the group
    for index, row in group.iterrows():
        new_time = (pd.to_datetime(row['real_time'], format='%H:%M:%S') + timedelta(minutes=offset)).time()
        group.at[index, 'real_time'] = new_time.strftime('%H:%M:%S')
    return group

# Apply the offset function to each group
data = data.groupby(data['date'].dt.date).apply(offset_time)

# Combine date and real_time to create datetime field
data['datetime'] = pd.to_datetime(data['date'].dt.date.astype(str) + ' ' + data['real_time'].astype(str))

# Sort the data by datetime
data = data.sort_values(by='datetime')

# Now that the data is sorted by datetime, reallocate sighting_ID as sequential numbers
data['sighting_ID'] = np.arange(1, len(data) + 1)

# Continue to save the simulated dataset as before
simulated_data_path = 'simulated-dugong-sightings-nGBR2023.csv'
data.to_csv(simulated_data_path, index=False)

print(f"Simulated dataset created and saved to {simulated_data_path}")
