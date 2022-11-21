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


    def position(self, z):
        # print(z)
        nearest_z = z[(np.abs(z - self.z_0)).argmin()]

        z_diff = nearest_z - self.z_0  # distance between starting point and first sensor
        multiplier = self.p_z/z_diff
        z_mult = multiplier * self.p_z
        x_mult = (multiplier * self.p_x) / z_mult # calculating the in momentum per z
        y_mult = (multiplier * self.p_y) / z_mult
        z_mult = z_mult / z_mult  # normalising z

        x_arr = z * x_mult
        y_arr = z * y_mult 

        # print(x_arr, y_arr)
        return x_arr, y_arr


        


if __name__ == "__main__":
    pass