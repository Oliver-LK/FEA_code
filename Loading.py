#***********************************************************
#  AUTHOR        : Oliver Clements
#  CREATE DATE   : 11/8/2023
#  PURPOSE       : ENME302 FEA code for distributed loads
#                  and point loads within a frame element
#  
# ***********************************************************


# Library imports
import numpy as np
from useful_funcs import remove
from Frames import Frame


# Transverse loading
def UDL_nodal_equivalent(frame: Frame, load: float):
    """ Gives the equivalent nodal load for a UDL"""

    nodal_equivalent = np.array([0, load * frame.length / 2, (load * frame.length ** 2 / 12), 0, load * frame.length / 2, (-load * frame.length ** 2 / 12)])
    # print(nodal_equivalent.T)
    return nodal_equivalent.T


def LVL_nodal_equivalent(frame: Frame, load: float):
    """ Gives the equivalent nodal load for a LVL"""

    nodal_equivalent = np.array([0, 3 * load * frame.length / 20, (load * frame.length ** 2 / 30), 0, 7 * load * frame.length / 20, (-load * frame.length ** 2 / 20)])
    return nodal_equivalent.T


def MSL_nodal_equivalent(frame: Frame, load: float):
    """ Gives the equivalent nodal load for a mid span load"""

    nodal_equivalent = np.array([0, load / 2, load * frame.length / 8, 0, load / 2, -load * frame.length / 8])
    return nodal_equivalent.T


def general_MSL_nodal_equivalent(frame: Frame, load: float, x: float):
    """ GIves the equivalent nodal forces for a load some distance along it"""

    nodal_equivalent = np.array([0, 
                                load * (1 - 3 * ((x ** 2) / (frame.length ** 2)) + 2 * ((x ** 3) / (frame.length ** 3))),
                                load * (x ** 3 / (frame.length ** 2) -  2 * (x ** 2 / frame.length) + x), 
                                0, 
                                load * (3 * ((x ** 2) / (frame.length ** 2)) - 2 * ((x ** 3) / (frame.length ** 3))), 
                                load * (x ** 3 / (frame.length ** 2) -  (x ** 2 / frame.length))])
    
    return nodal_equivalent.T



# Axial loading
def DAL_nodal_equivalent(frame: Frame, load: float):
    """ Gives the equivalent nodal load for a uniform distributed axial load"""

    nodal_equivalent = np.array([load * frame.length / 2, 0, 0, load * frame.length / 2, 0, 0])
    return nodal_equivalent.T


def CAL_nodal_equivalent(frame: Frame, load: float, distance: float):
    """ Gives the equivalent nodal load an axial load somewhere along the frame element"""

    nodal_equivalent = load * np.array([1 - (distance / frame.length), 0, 0, distance / frame.length, 0, 0])
    return nodal_equivalent

# Transforms
def load_to_global(frame: Frame, element_force: np.ndarray):
    """ Converts the load into global coordinates"""

    up_lambda = frame.give_transform_matrix()
    global_force = up_lambda.T @ element_force

    return global_force


def give_equivalent_force(frame: Frame, global_force: np.ndarray):
    """ Gives the equivalent nodal force from a distributed load in global coords"""

    global_coord_force = frame.assembly_matrix @ global_force
    return global_coord_force


def do_DL_steps(frame: Frame, element_force: float):
    """ Does the steps for you"""

    global_force = load_to_global(frame, element_force)
    global_coord_force = give_equivalent_force(frame, global_force)
    return global_coord_force
