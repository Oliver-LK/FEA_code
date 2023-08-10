import numpy as np

MAGNITUDE_FACTOR = 1 * 10 ** 6

def give_area_of_circle(diameter):
    return (np.pi * diameter ** 2) / 4

def give_mm_to_m(measurement):
    """ Converts mm to m"""
    return measurement / 1000

def remove(matrix: np.ndarray):
    """ Determines if array is 1d or 2d and calls the appropriate method"""
    matrix_shape = matrix.shape

    if len(matrix_shape) == 2:
        matrix = remove_insignificant_values_2d(matrix)

    elif len(matrix_shape) == 1:
        matrix = remove_insignificant_values_1d(matrix)

    return matrix

def remove_insignificant_values_1d(matrix: np.ndarray):
    """ Removes values if from a 2D array if they are insignificant"""
    row_index = 0
    
    largest_value = np.max(matrix)

    while row_index < len(matrix):

        if(abs(matrix[row_index]) < (largest_value / MAGNITUDE_FACTOR)):
            matrix[row_index] = 0

        row_index += 1
    
    return matrix

def remove_insignificant_values_2d(matrix: np.ndarray):
    """ Removes values if from a 2D array if they are insignificant"""
    row_index = 0
    col_index = 0
    
    largest_value = np.max(matrix)

    while col_index < len(matrix):
        while row_index < len(matrix[0]):

            if((abs(matrix[col_index][row_index]) < (largest_value / MAGNITUDE_FACTOR) or (abs(matrix[col_index][row_index] < 1 * 10 ** -15)))):
                matrix[col_index][row_index] = 0

            row_index += 1

        col_index += 1
        row_index = 0
    
    return matrix

