import numpy as np
import matplotlib.pyplot as plt
from useful_funcs import *
from Bars import *

# Assembly matrices
ASSEM1 = np.array([[0, 0, 1, 0], [0, 0, 0, 1]])
ASSEM2 = np.array([[0, 0, 1, 0], [0, 0, 0, 1]])

# External force
Q = np.array([0, 100 * 10 ** 3])

# Define bars
bar_1 = Bar(10, give_area_of_circle(0.1), 200, 0, ASSEM1)
bar_2 = Bar(10 * np.sqrt(2), give_area_of_circle(0.1), 200, 45, ASSEM2)

# Bars in list
bar_list = [bar_1, bar_2]

# Global stiffness of each local element
k1G = bar_1.give_global_contribution_stiffness()
k2G = bar_2.give_global_contribution_stiffness()

# Total global stiffness matrix
total_k = k1G + k2G

# Global deflections
global_deflections = give_global_deflections(total_k, Q)

# Forces in each element
forces1 = bar_1.give_nodal_forces(global_deflections)
forces2 = bar_2.give_nodal_forces(global_deflections)

# Print all of interest
print(forces1)
print(forces2)
print(bar_1)

