'''Creating a particle class'''


import numpy as np


class particle:
    def __init__(self, momentum, start):
        '''Constructor for particle class
        
        Args:
        momentum: tuple/list - a momentum vector of the form (p_x, p_y, p_z)
        start: tuple/list - the starting position of the particle that we measure (x, y, z))
        '''

        # print(momentum)
        self.p_arr = np.array(momentum)
        self.p_x, self.p_y, self.p_z = momentum
        self.x_0, self.y_0, self.z_0 = start


    def calibrate(self, z_arr):
        # print(z)
        nearest_z = z_arr[(np.abs(z_arr - self.z_0)).argmin()]

        z_diff = nearest_z - self.z_0  # distance between starting point and first sensor
        multiplier = self.p_z/z_diff
        z_mult = multiplier * self.p_z
        self.x_mult = (multiplier * self.p_x) / z_mult # calculating the in momentum per z
        self.y_mult = (multiplier * self.p_y) / z_mult
        self.z_mult = z_mult / z_mult  # normalising z

        self.z_arr = z_arr
        self.x_arr = z_arr * self.x_mult
        self.y_arr = z_arr * self.y_mult 

        # print(x_arr, y_arr)


    def position(self, z_value):
        # print(z_value)
        index = np.where(z_value == self.z_arr)
        # print(self.x_arr[index], self.y_arr[index], self.z_arr[index])
        return self.x_arr[index], self.y_arr[index], self.z_arr[index] 

        


if __name__ == "__main__":
    pass