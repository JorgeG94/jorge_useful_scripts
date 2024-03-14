import pandas as pd
import matplotlib.pyplot as plt

# Read the data from the file
file_path = '/home/jorgegv/work/scripts_jorge/get_ordered_dimer/distances.txt'  # Replace with the actual path to your file
data = pd.read_table(file_path, skiprows=2, delim_whitespace=True, header=None)
data.columns = ['ID', 'E', 'dE', 'dimer_distance']

# Sort the DataFrame based on the 'dimer_distance' column
data_sorted = data.sort_values(by='dimer_distance')


# Plot dE against dimer_distance
plt.plot(data_sorted['dimer_distance'], data_sorted['dE'], marker='o', linestyle='-')
plt.xlabel('Dimer Distance')
plt.ylabel('dE')
plt.title('Plot of dE against Dimer Distance')

# Identify points where abs(dE) < 1E-6
threshold = 3.8E-7
filtered_data = data_sorted[abs(data_sorted['dE']) < threshold]

for distance in filtered_data['dimer_distance']:
    print("distance ", distance)

# Plot lines at points where abs(dE) < 1E-6
for distance in filtered_data['dimer_distance']:
    plt.axvline(x=distance, color='r', linestyle='--')
    plt.text(distance, 0, f'{distance:.2f}', rotation=90, va='bottom', ha='right', color='r')

# Save the plot to a PNG file
plt.savefig('output_plot.png')