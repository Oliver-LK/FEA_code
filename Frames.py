#***********************************************************
#  AUTHOR        : Oliver Clements
#  CREATE DATE   : 2/8/2023
#  PURPOSE       : ENME302 FEA code for FRAME ELEMENTS. This can solve
#                  the global and nodal deflections, forces a
#                  simple system. It appears it only allows
#                  loads to be placed at point of freedom in
#                  structure
#  
# ***********************************************************


# Library imports
import numpy as np
import matplotlib.pyplot as plt
from useful_funcs import remove

N_POINTS = 25       # Num of points in each element
PLOT_MAGNIFY = 30  # Magnification factor to actually see the displacements


class Frame:

    def __init__(self, length: float, area: float, second_MoA: float, E: float, angle: float, assembly_matrix: np.ndarray, element_nodal_deflections: np.ndarray, global_nodal_deflections: np.ndarray):
        """ Defines the parameters of a frame element
            length: m
            second_MoA: m ** 4
            area: m^2
            E: GPa
            angle is angle from global frame. Defined in degrees
            assembly matrix: Matrix that defines where the deflections occur within the structure adn relates this to the individual element
            element_nodal_deflections: Element deflections in element coordinates
            
        """
        
        self.length = length
        self.area = area
        self.second_MoA = second_MoA
        self.E = E * 10 ** 9
        self.angle = angle
        self.assembly_matrix = assembly_matrix
        self.element_nodal_deflections = element_nodal_deflections
        self.global_nodal_deflections = global_nodal_deflections

        self.beta = (area * length ** 2) / second_MoA


    def __str__(self):
        return f"Length: {self.length}m, Area {self.area} m^2, Second moment of area {self.second_MoA} m^4 E: {self.E} GPa\n{self.angle} degrees from global reference frame and has the following assembly matrix \n{self.assembly_matrix}"
    

    def give_local_stiffness(self):
        """ Gives the locall stiffness matrix of element"""
        
        matrix = np.array([[self.beta, 0, 0, -self.beta, 0, 0],
                           [0, 12, 6 * self.length, 0, -12, 6 * self.length],
                           [0, 6 * self.length, 4 * self.length ** 2, 0, -6 * self.length, 2 * self.length ** 2],
                           [-self.beta, 0, 0, self.beta, 0, 0],
                           [0, -12, -6 * self.length, 0, 12, -6 * self.length],
                           [0, 6 * self.length, 2 * self.length ** 2, 0, -6 * self.length, 4 * self.length ** 2]])
        
        return matrix * (self.second_MoA * self.E / (self.length ** 3))
    

    def give_transform_matrix(self):
        """ Gives the horrid transform matrix"""
        angle = np.deg2rad(self.angle)

        # Building the matrix
        zero_matrix = np.zeros((3, 3))
        non_zero_part = np.array([[np.cos(angle), np.sin(angle), 0],
                                  [-1 * np.sin(angle), np.cos(angle), 0],
                                  [0, 0, 1]])
        
        up_lambda = np.concatenate((np.concatenate((non_zero_part, zero_matrix), axis=1),np.concatenate((zero_matrix, non_zero_part), axis=1)), axis=0)

        return up_lambda
        

    def give_global_stiffness(self):
        """ Gives the global stiffness of the element"""
        
        up_lambda = Frame.give_transform_matrix(self)
        # Maths part
        local_stiffness = Frame.give_local_stiffness(self)

        global_stiffness = up_lambda.T @ local_stiffness @ up_lambda

        return (global_stiffness) # Could be a source of error
    

    def give_global_contribution_stiffness(self):
        """ Gives the element stiffness contribution to global frame
            assembly_matrix: found manually """
        
        return self.assembly_matrix @ Frame.give_global_stiffness(self) @ self.assembly_matrix.T
    
    
    def give_element_nodal_forces(self, element_deflections: np.ndarray):
        """ Gives the nodal forces in element coordinates"""
        local_stiffness = Frame.give_local_stiffness(self)

        return local_stiffness @ element_deflections


    def give_global_nodal_forces(self, global_deflections: np.ndarray):
        """ Gives the nodal forces in global coordinates"""
        global_stiffness = Frame.give_global_stiffness(self)

        return global_stiffness @ (self.assembly_matrix.T @ global_deflections)
    

    def give_nodal_displacements_global_coord(self, global_deflections: np.ndarray):
        """Gives the displacements of each node in global_coordinates"""
        deflection = self.assembly_matrix.T @ global_deflections
        self.global_nodal_deflections = deflection
    

    def give_nodal_displacements_element_coord(self, global_deflections: np.ndarray):
        """ Gives the displacements of each node in element coordinates"""
        
        up_lambda = Frame.give_transform_matrix(self)
        
        displacements = up_lambda @ self.assembly_matrix.T @ global_deflections
        self.element_nodal_deflections = displacements



def give_overall_stiffness_matrix(element_list: list[Frame], num_elements):
    """ Returns the total global stiffness matrix of all elements"""

    total_global_stiffness = np.zeros((num_elements + 1, num_elements + 1))

    for element in element_list:
        total_global_stiffness += element.give_global_contribution_stiffness()

    return total_global_stiffness


def give_global_deflections(overall_global_stiffness: np.ndarray, force: np.ndarray):
    """ Finds the global deflections"""        

    deflections = np.linalg.inv(overall_global_stiffness) @ force
    return deflections
    


class Plot_frame(Frame):
    
    def __init__(self, length: float, area: float, second_MoA: float, E: float, angle: float, assembly_matrix: np.ndarray, element_nodal_deflections: np.ndarray, global_nodal_deflections: np.ndarray, global_coord: dict, axial_shape: np.ndarray, transverse_shape: np.ndarray, element_axial_d, element_transverse_d, xg_deflection, yg_deflection, xg_baseline, yg_baseline):
        """ Defines the necessary plotting elements for the frame element
            global_coord: (x, y) in m

            All these are defined later:
            axial_shape is a tuple which is (phi1, phi2)
            transverse_shape: N is an array like (n1, n2, 3, n4)
            element_axial_d: Axial deflection
            element_transverse_d: Transverse deflection


        """
        super().__init__(length, area, second_MoA, E, angle, assembly_matrix, element_nodal_deflections, global_nodal_deflections)

        self.global_coord = global_coord
        self.axial_shape = axial_shape
        self.transverse_shape = transverse_shape
        self.element_axial_d = element_axial_d
        self.element_transverse_d = element_transverse_d
        self.xg_deflection = xg_deflection
        self.yg_deflection = yg_deflection
        self.xg_baseline = xg_baseline
        self.yg_baseline = yg_baseline


    def __str__(self):
        return f"Length: {self.length}m, Area {self.area} m^2, Second moment of area {self.second_MoA} m^4 E: {self.E} GPa\n{self.angle} degrees from global reference frame and has the following assembly matrix \n{self.assembly_matrix}"


    def find_shape_functions(self):
        """ Finds the axial and transverse shape functions and assigns it to the class"""
        x_arr = np.linspace(0, self.length, N_POINTS)

        # Defining axial shape functions
        phi1, phi2 = np.zeros((N_POINTS)), np.zeros((N_POINTS))

        # Defining transverse shape functions
        n1, n2, n3, n4 = np.zeros((N_POINTS)), np.zeros((N_POINTS)), np.zeros((N_POINTS)), np.zeros((N_POINTS))

        for i, x in enumerate(x_arr):
            phi1[i] = 1 - (x / self.length)
            phi2[i] = x / self.length

            n1[i] = 1 - ((3 * x ** 2) / self.length ** 2) + ((2 * x ** 3) / self.length ** 3)
            n2[i] = ((x ** 3) / self.length ** 2) - ((2 * x ** 2) / self.length) + x
            n3[i] = ((3 * x ** 2) / self.length ** 2) - ((2 * x ** 3) / self.length ** 3)
            n4[i] = ((x ** 3) / self.length ** 2) - ((x ** 2) / self.length)

        self.axial_shape = (phi1, phi2)
        self.transverse_shape = np.column_stack((n1, n2, n3, n4))


    def find_axial_displacements(self, d1: float, d4: float):
        """ Finds the axial_displacements and assigns it to the class"""
        axial_displacements = (self.axial_shape[0] * d1) + (self.axial_shape[1] * d4)
        self.element_axial_d = axial_displacements

    
    def find_transverse_displacements(self, displacements: np.ndarray):
        """ Finds the transverse displacements and assigns it to the class"""

        transverse_displacements = self.transverse_shape @ displacements.T
        self.element_transverse_d = transverse_displacements

    
    def find_baselines(self):
        """ Finds the baselines of each element. This is essentially the undflected shape"""

        self.xg_baseline = np.linspace(self.global_coord["X1G"], self.global_coord["X2G"], N_POINTS)
        self.yg_baseline = np.linspace(self.global_coord["Y1G"], self.global_coord["Y2G"], N_POINTS)


    def find_deflections_XG_YG(self):
        """ Finds the deflections in x of a frame"""
        xg_deflection = np.zeros((N_POINTS))
        yg_deflection = np.zeros((N_POINTS))

        angle = np.deg2rad(self.angle)

        for point in range(0, len(self.element_axial_d)):
            xg_deflection[point] = ((self.element_axial_d[point] * np.cos(angle)) - (self.element_transverse_d[point] * np.sin(angle))) * PLOT_MAGNIFY + self.xg_baseline[point]
            yg_deflection[point] = ((self.element_axial_d[point] * np.sin(angle)) + (self.element_transverse_d[point] * np.cos(angle))) * PLOT_MAGNIFY + self.yg_baseline[point]

        self.xg_deflection = xg_deflection
        self.yg_deflection = yg_deflection



def copy_class_details(class_copy: Frame):
    """ Copy's details from Frame class to a new Plot_frame class"""
    new_class = Plot_frame(length=class_copy.length, area=class_copy.area, second_MoA=class_copy.second_MoA, angle=class_copy.angle, E=class_copy.angle, assembly_matrix=class_copy.assembly_matrix, element_nodal_deflections=class_copy.element_nodal_deflections, global_nodal_deflections= class_copy.global_nodal_deflections, global_coord=None, axial_shape=None, transverse_shape=None,element_axial_d=None, element_transverse_d=None, xg_deflection = None, yg_deflection = None, xg_baseline = None, yg_baseline = None)
    
    return new_class


def plot_deflection(frame_list: list[Plot_frame]):
    """ Plots the undeflected and the deflected frame"""

    axes = plt.axes()   
    axes.grid(True)
    axes.set_title(f"Deflected vs Undeflected structure using\n a magnification of {PLOT_MAGNIFY}")
    axes.set_xlabel("X (m)")
    axes.set_ylabel("Y (m)")   
    
    for frame in frame_list:
        axes.plot(frame.xg_baseline, frame.yg_baseline,'b-')
        axes.plot(frame.xg_deflection, frame.yg_deflection,'r-')
        
    axes.legend(["Undeflected Frame", "Deflected Frame"])
    plt.show()         
