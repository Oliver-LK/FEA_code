import numpy as np
import matplotlib.pyplot as plt
from useful_funcs import *
from Frames import *


# Necessary lines
frame1 = Frame(length=10, area=1 * 10 ** -5, second_MoA=5 * 10 ** -4, angle=-90,  E=200, assembly_matrix=None, element_nodal_deflections=None, global_nodal_deflections=None)
frame2 = Frame(length=10, area=1 * 10 ** -5, second_MoA=5 * 10 ** -4, angle=0,  E=200, assembly_matrix=None, element_nodal_deflections=None, global_nodal_deflections=None)

frame_list = [frame1, frame2]

frame1.assembly_matrix = np.array([[0, 0, 0, 1, 0, 0],
                                   [0, 0, 0, 0, 0, 1],
                                   [0, 0, 0, 0, 0, 0]])

frame2.assembly_matrix = np.array([[1, 0, 0, 0, 0, 0],
                                   [0, 0, 1, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 1]])

# Non Necessary
local_stiffness1 = frame1.give_local_stiffness()
global_stiffness1 = frame2.give_global_contribution_stiffness()


# Necessary lines
overall_global_stiffness = give_overall_stiffness_matrix(frame_list)

# External forcing term
Q = np.array([0, 140000, 0])
Q = Q.T         # Transpose so its a vector

overall_deflections = give_global_deflections(overall_global_stiffness, Q)
# print(overall_deflections)

# Element 1
frame1.give_nodal_displacements_global_coord(overall_deflections)
frame1.give_nodal_displacements_element_coord(overall_deflections)
frame1.give_global_nodal_forces(overall_deflections)
element_nodal_forces1 = frame1.give_element_nodal_forces(frame1.element_nodal_deflections)

# Element 2
frame2.give_nodal_displacements_global_coord(overall_deflections)
frame2.give_nodal_displacements_element_coord(overall_deflections)
frame2.give_global_nodal_forces(overall_deflections)
element_nodal_forces2 = frame2.give_element_nodal_forces(frame2.element_nodal_deflections)

print(remove(frame1.element_nodal_deflections))
print(remove(frame2.element_nodal_deflections))



# Plotting stuff
# Need
# Copying class details for new child class
frame_plot_list = []

for frame in frame_list:
    frame_plot_list.append(copy_class_details(frame))


# Set global coordinates of each node of each element
frame_plot_list[0].global_coord = {"X1G": 0, 'X2G': 0, "Y1G": frame1.length, "Y2G": 0}
frame_plot_list[1].global_coord = {"X1G": 0, 'X2G': frame2.length, "Y1G": 0, "Y2G": 0}

for frame in frame_plot_list:
    frame.find_shape_functions()
    element_nodal_deflections1 = frame.give_nodal_displacements_element_coord(overall_deflections)
    frame.find_baselines()
    frame.find_axial_displacements(frame.element_nodal_deflections[0], frame.element_nodal_deflections[3])
    frame.find_transverse_displacements(np.array([frame.element_nodal_deflections[1], frame.element_nodal_deflections[2], frame.element_nodal_deflections[4], frame.element_nodal_deflections[5]]))
    frame.find_deflections_XG_YG()
    print(f"Axial deformation: {remove(frame.element_axial_d)}")
    print(f"Transverse deformation: {remove(frame.element_transverse_d)}\n")



# Plot
plot_deflection(frame_plot_list)

# frame1_plot.find_shape_functions()
# frame2_plot.find_shape_functions()
# # print(frame1_plot.transverse_shape)

# # print(remove_insignificant_values(element_nodal_deflections1))
# element_nodal_deflections1 = frame1.give_nodal_displacements_element_coord(overall_deflections)
# element_nodal_deflections2 = frame2.give_nodal_displacements_element_coord(overall_deflections)


# # Need this code
# frame1_plot.find_baselines()
# frame2_plot.find_baselines()

# frame1_plot.find_axial_displacements(np.array([element_nodal_deflections1[0]]), np.array([element_nodal_deflections1[3]]))
# frame2_plot.find_axial_displacements(np.array([element_nodal_deflections2[0]]), np.array([element_nodal_deflections2[3]]))

# frame1_plot.find_transverse_displacements(np.array([element_nodal_deflections1[1], element_nodal_deflections1[2], element_nodal_deflections1[4], element_nodal_deflections1[5]]))
# frame2_plot.find_transverse_displacements(np.array([element_nodal_deflections2[1], element_nodal_deflections2[2], element_nodal_deflections2[4], element_nodal_deflections2[5]]))

# frame1_plot.find_deflections_XG_YG()
# frame2_plot.find_deflections_XG_YG()



# print(frame1_plot.element_axial_d)
# print(frame2_plot.element_axial_d)


# print(frame1_plot.element_transverse_d)
# print(frame2_plot.element_transverse_d)

# print(len(frame1_plot.xg_baseline))
# print(len(frame1_plot.xg_deflection))
