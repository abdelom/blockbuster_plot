import matplotlib.pyplot as plt
import numpy as np
import math

# Function to read y values from the file
def read_y_values_from_file(file_path):
    with open(file_path, 'r') as file:
        # Read lines and convert values to a list of floats
        y_values = [list(map(float, line.strip().split())) for line in file.readlines()]
    return y_values

# Function to calculate consecutive differences
def calculate_consecutive_differences(y_values):
    return [[y[i + 1] - y[i] for i in range(len(y) - 1)] for y in y_values]

# Function to calculate the sum of differences for each x value
def calculate_sum_of_differences(differences):
    max_length = max(len(diff) for diff in differences)  # Handle different lengths
    summed_diff_values = np.array([sum(diff[i] for diff in differences if i < len(diff))  for i in range(max_length)])
    return summed_diff_values

# Function to calculate cumulative sum of the sums
def calculate_cumulative_sum(summed_differences):
    return np.cumsum(summed_differences)

# Path to the file (replace with your file path)
file_path = 'grid.txt'

# Read y values from the file
y_values = read_y_values_from_file(file_path)

# Calculate consecutive differences for each line
differences = calculate_consecutive_differences(y_values)

# Calculate the sum of differences for each x value
summed_differences = calculate_sum_of_differences(differences)

# Calculate the cumulative sum of the summed differences
cumulative_sum = calculate_cumulative_sum(summed_differences)

# x values (assuming the same range for each curve)
x = []
i = 0
while i < len(y_values[0]):
    x.append((i+1) * 4e-2)
    i += 1

print(x, len(x))
x_diff = []
i = 0
while i < len(y_values[0]) - 1:
     x_diff.append((i+1) * 4e-2)
     i += 1
#x = range(1, len(y_values[0]) + 1)
#x_diff = range(1, len(y_values[0]))  # For differences, one less value
x_cumsum = range(1, len(cumulative_sum) + 1)  # For cumulative sums

# Create the first plot (original curves)
# Compute differences for variation plot
differences = [np.diff(y) for y in y_values]

plt.figure(figsize=(10, 8))
for idx, y in enumerate(y_values):
    if idx > 5:
        break
    plt.plot(x, y, marker='o', label=f'Length of branches supporting {idx + 1} leaves')
plt.xlabel('Time in Ne generations')
plt.ylabel('Branch length in Ne generations')
plt.title('Cumulative length of branches supporting a certain number of leaves over time')
plt.legend()
plt.savefig('cumulative_branch_lengths.png')  # Save the first plot
plt.show()

# Plot 2: Branch Length Differences
plt.figure(figsize=(10, 8))
x_diff = x[:-1]  # Exclude the last point for differences
for idx, diff in enumerate(differences):
    if idx > 5:
        break
    plt.plot(x_diff, diff, marker='x', label=f'Variation for branches supporting {idx + 1} leaves')
plt.xlabel('Time in Ne generations')
plt.ylabel('Branch length in Ne generations')
plt.title('Variation in branch lengths supporting a certain number of leaves over time')
plt.legend()
plt.savefig('branch_length_differences.png')  # Save the second plot
plt.show()

# Show legend


# Create the third plot for the sum of differences
plt.figure(figsize=(10, 8))

def f(x):
    return np.log( np.exp(-x) / x)
# Plot the summed differences
plt.plot(range(1, len(summed_differences)), summed_differences[:-1], marker='x',color='purple', label='Somme des différences')
plt.ylim(0.0, 0.5)
#plt.plot(range(1, len(summed_differences)), f(np.arange(0.1, 10.1, 0.1)), label=r'$f(x) = \frac{e^{-x}}{x}$')
# Add labels and title
plt.xlabel('X values')
plt.ylabel('Somme des différences')
plt.title('Somme des différences pour chaque valeur de X')

# Show legend
plt.legend()

# Show the summed differences plot
plt.show()

# Create the fourth plot for the cumulative sum
plt.figure(figsize=(10, 8))

# Plot the cumulative sum of summed differences
plt.plot(x_cumsum, cumulative_sum, color='orange', marker='x', label='Somme cumulée des différences')
# Add labels and title
plt.xlabel('X values')
plt.ylabel('Somme cumulée des différences')
plt.title('Somme cumulée des différences pour chaque valeur de X')

# Show legend
plt.legend()

# Show the cumulative sum plot
plt.show()
