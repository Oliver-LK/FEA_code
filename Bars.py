#***********************************************************
#  AUTHOR        : Oliver Clements
#  CREATE DATE   : 31/7/2023
#  PURPOSE       : ENME302 FEA code for BAR ELEMENTS. This can solve
#                  the global and nodal deflections, forces a
#                  simple system. It appears it only allows
#                  loads to be placed at point of freedom in
#                  structure
#  
#  ***********************************************************

# Library imports
import numpy as np
import matplotlib.pyplot as plt


class Bar:
    """ Class for bar with associated methods"""

    def __init__(self, length: float, area: float, E: float, angle: float, assembly_matrix: np.ndarray):
        """ Defines the parameters of a bar
            length: m
            area: m^2
            E: GPa
            angle is angle from global frame
        """
        
        self.length = length
        self.area = area
        self.E = E * 10 ** 9
        self.angle = angle
        self.assembly_matrix = assembly_matrix

    def __str__(self):
        return f"Length: {self.length}m, Area {self.area} m^2, E: {self.E} GPa\n{self.angle} degrees from global reference frame and has the following assembly matrix \n{self.assembly_matrix}"
    
    # Methods
  
    def give_local_stiffness(self):
        """ Gives the stiffness matrix"""
        return (self.area * self.E /self.length) * np.array([[1, -1],[-1, 1]])
    
    
    def give_global_stiffness(self):
        """ Gives the stiffness in global terms
            angle: degrees"""
        angle = np.deg2rad(self.angle)

        up_lambda = np.array([[np.cos(angle), np.sin(angle), 0, 0],
                             [0, 0, np.cos(angle), np.sin(angle)]])
        
        local_stiffness = Bar.give_local_stiffness(self)

        return up_lambda.T @ local_stiffness @ up_lambda
    
    
    def give_global_contribution_stiffness(self):
        """ Gives the element stiffness contribution to global frame
            assembly_matrix: found manually """
        
        return self.assembly_matrix @ Bar.give_global_stiffness(self) @ self.assembly_matrix.T
    
    
    def give_global_nodal_forces(self, global_deflections: np.ndarray):
        """ Gives the nodal forces"""
        local_stiffness = Bar.give_global_stiffness(self)

        return local_stiffness @ (self.assembly_matrix.T @ global_deflections)
    
    
    def give_nodal_displacements(self, global_deflections: np.ndarray):
        """ Gives the displacements of each node"""
        angle = np.deg2rad(self.angle)

        up_lambda = np.array([[np.cos(angle), np.sin(angle), 0, 0],
                             [0, 0, np.cos(angle), np.sin(angle)]])
        
        displacements = up_lambda @ self.assembly_matrix.T @ global_deflections

        return displacements
    
    
    def give_element_strain(self, global_deflections: np.ndarray):
        """ Gives the strain of the element"""
        displacements = Bar.give_nodal_displacements(self, global_deflections)
        strain = (displacements[1] - displacements[0]) / self.length

        return strain
    
    
    def give_nodal_forces(self, global_deflections: np.ndarray):
        """ Gives the forces at each node of an element. This is a vector with 2 elements"""

        local_stiffness = Bar.give_local_stiffness(self)
        nodal_displacements = Bar.give_nodal_displacements(self, global_deflections)
        forces = local_stiffness @ nodal_displacements

        return forces


def give_overall_stiffness_matrix(element_list: list[Bar]):
    """ Returns the total global stiffness matrix of all elements"""

    total_global_stiffness = np.array([])

    for element in element_list:
        total_global_stiffness += Bar.element.give_global_contribution_stiffness()

    return total_global_stiffness


def give_global_deflections(global_stiffness: np.ndarray, force:np.ndarray):
    """ Finds the global deflections
    """        

    deflections = np.linalg.inv(global_stiffness) @ force
    return deflections
    
