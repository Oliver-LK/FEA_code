# Example code for distributed loads on frame elements

import numpy as np
from useful_funcs import *
from Frames import *
from Loading import *

frame1 = Frame(length=3, area=2 * 10 ** -4, second_MoA=6 * 10 ** -6, angle=80,  E=30, assembly_matrix=None, element_nodal_deflections=None, global_nodal_deflections = None)
frame2 = Frame(length=2, area=2 * 10 ** -4, second_MoA=6 * 10 ** -6, angle=0,  E=30, assembly_matrix=None, element_nodal_deflections=None, global_nodal_deflections = None)


frame_list = [frame1, frame2]

frame1.assembly_matrix = np.array([[0, 0, 0, 1, 0, 0],
                                   [0, 0, 0, 0, 1, 0],
                                   [0, 0, 0, 0, 0, 1],
                                   [0, 0, 0, 0, 0, 0]])

frame2.assembly_matrix = np.array([[1, 0, 0, 0, 0, 0],
                                   [0, 1, 0, 0, 0, 0],
                                   [0, 0, 1, 0, 0, 0],
                                   [0, 0, 0, 0, 1, 0]])



# Non Necessary
local_stiffness1 = frame1.give_local_stiffness()
global_stiffness1 = frame1.give_global_contribution_stiffness()
global_stiffness2 = frame2.give_global_contribution_stiffness()

print(f"Global stiffness for element 1:\n{remove(global_stiffness1)}\n")
print(f"Global stiffness for element 2:\n{remove(global_stiffness2)}\n")



# Necessary lines
overall_global_stiffness = give_overall_stiffness_matrix(frame_list, frame1.assembly_matrix.shape[0])
print(f"Global_stiffness matrix:\n{np.around(overall_global_stiffness, decimals=0)}\n")


# External UDL
f_eq_udl = UDL_nodal_equivalent(frame2, -10000)
print(f"UDL:\n{f_eq_udl}")
Q_eq_udl = do_DL_steps(frame2, f_eq_udl)
Q_eq_udl.T


# Total force
Q =  Q_eq_udl
print(f"Q:\n{Q}")

overall_deflections = give_global_deflections(overall_global_stiffness, Q)
print(f"Global Deflections:\n{np.around(overall_deflections, decimals=5)}\n")

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
# print(f"Element nodal forces are:\n{element_nodal_forces2}")



print(f"Support force for 1:\n{element_nodal_forces1}\n")
print(f"Support force for 2:\n{element_nodal_forces2}\n")


axial_stress = frame2.give_axial_stress(element_nodal_forces2[0])
print(f"Axial stress:\n{axial_stress}")

strain1 = frame1.give_strain()
print(f"Element strain: {strain1}")

# Plotting stuff
# Copying class details for new child class
frame_plot_list = []

for frame in frame_list:
    frame_plot_list.append(copy_class_details(frame))


# Set global coordinates of each node of each element
frame_plot_list[0].global_coord = {"X1G": 0, 'X2G': frame1.length * np.cos(np.deg2rad(frame1.angle)), "Y1G": 0, "Y2G": frame1.length * np.sin(np.deg2rad(frame1.angle))}
frame_plot_list[1].global_coord = {"X1G": frame1.length * np.cos(np.deg2rad(frame1.angle)), 'X2G': frame2.length + frame1.length * np.cos(np.deg2rad(frame1.angle)), "Y1G": frame1.length * np.sin(np.deg2rad(frame1.angle)), "Y2G": frame1.length * np.sin(np.deg2rad(frame1.angle))}
for frame in frame_plot_list:
    frame.find_shape_functions()
    element_nodal_deflections1 = frame.give_nodal_displacements_element_coord(overall_deflections)
    frame.find_baselines()
    frame.find_axial_displacements(frame.element_nodal_deflections[0], frame.element_nodal_deflections[3])
    frame.find_transverse_displacements(np.array([frame.element_nodal_deflections[1], frame.element_nodal_deflections[2], frame.element_nodal_deflections[4], frame.element_nodal_deflections[5]]))
    frame.find_deflections_XG_YG()
    # print(f"Axial deformation: {remove(frame.element_axial_d)}")
    # print(f"Transverse deformation: {remove(frame.element_transverse_d)}\n")

# Plot
plot_deflection(frame_plot_list)
