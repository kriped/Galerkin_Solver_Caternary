import numpy as np

import matplotlib.pyplot as plt

# Load data from results.out
data = np.loadtxt('results.out')

# Extract X and Y components
Y = data  # Assuming the first column is X
X = np.linspace(-1, 1, len(Y))

# Create the plot
plt.figure()
plt.plot(X, Y, marker='o', linestyle='-', color='b', label='Vector Data')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Plot of Calculation Data')
plt.legend()
plt.grid(True)

# Show the plot
plt.show()