'''Crreating classes for the left and right side VELO sensors'''


import matplotlib.pyplot as plt
import numpy as np
from itertools import product, combinations



# We have 26 sensors for each side. 
z_left_sensors = np.array([-277, -252, -227, -202, -132, -62, -37, -12, 13, 38, 63, 88, 
113, 138, 163, 188, 213, 238, 263, 325, 402, 497, 616, 661, 706, 751])

z_left_sensors_bounds = np.array(list(zip(z_left_sensors, z_left_sensors + 0.2)))

z_right_sensors = np.array([-289, -264, -239, -214, -144, -74, -49, -24, 1, 26, 51, 76, 
101, 126, 151, 176, 201, 226, 251, 313, 390, 485, 604, 649, 694, 739])

def plot_sensor(x, y, axis):
    # Right Detector
    axis.plot([5.1, 33.26], [-5.1, -5.1], color="red")
    axis.plot([5.1, 33.26], [33.51, 33.51], color="red")
    axis.plot([5.1,5.1], [-5.1,33.51], color="red")
    axis.plot([33.26, 33.26], [-5.1, 33.51], color="red")
    axis.plot([-10.2, 37.41], [-33.51,-33.51], color="red")
    axis.plot([-5.1, 37.41], [-5.1,-5.1], color="red")
    axis.plot([-5.1, -5.1], [-33.51, -5.1], color="red")
    axis.plot([37.41, 37.41], [-33.51, -5.1], color="red")

    # Left Detector
    axis.plot([-5.1, -33.26], [5.1, 5.1], color="blue")
    axis.plot([-5.1, -33.26], [-33.51, -33.51], color="blue")
    axis.plot([-5.1, -5.1], [5.1, -33.51], color="blue")
    axis.plot([-33.26, -33.26], [5.1, -33.51], color="blue")
    axis.plot([5.1, -37.41], [33.51,33.51], color="blue")
    axis.plot([5.1, -37.41], [5.1,5.1], color="blue")
    axis.plot([5.1, 5.1], [33.51, 5.1], color="blue")
    axis.plot([-37.41, -37.41], [33.51, 5.1], color="blue")

    axis.scatter(x, y, label="particle", color="green")
    axis.legend(loc="upper right")
    axis.set_xlabel("x position (mm)")
    axis.set_ylabel("y postion (mm)")
    axis.set_title("Plot of particle incidence on sensor.")
    

class right_sensor():
    def __init__(self, z_pos):
        self.z_pos = z_pos


    def in_sensor(self, x, y):
        '''Calculating if the particle lies within the right hand sensor.
        
        Args: 
        x: float - x position of particle.
        y: float - y position of particle.
        
        Returns: 
        1 or 0: returns 1 if inside the sensor, or 0 if outside.'''

        index = np.where((((5.1 < x) & (x < 33.26)) & ((-5.1 < y) & (y < 33.51))) | (((-5.1 < x) & (x < 37.41)) & ((-33.51 < y) & (y < -5.1))))
        hits_y, hits_x = y[index], x[index]

        return hits_y, hits_x


class left_sensor():
    def __init__(self, z_pos):
        self.z_pos = z_pos


    def in_sensor(self, x, y):
        '''Calculating if the particle lies within the left hand sensor.
        
        Args: 
        x: float - x position of particle.
        y: float - y position of particle.
        
        Returns: 
        1 or 0: returns 1 if inside the sensor, or 0 if outside.'''
        # print(x)
        index = np.where((((-33.26 < x) & (x < -5.1)) & ((-33.51 < y) & (y < 5.1))) | (((-37.41 < x) & (x < 5.1)) & ((5.1 < y) & (y < 33.51))))
        hits_y, hits_x = y[index], x[index]
        
        return hits_y, hits_x
        

if __name__ == "__main__":
    print(right_sensor(-272).in_sensor(20, 67))  # should produce outside
    print(right_sensor(-272).in_sensor(-7, -20))  # should produce inside
    print(left_sensor(-272).in_sensor(-7, -20))
    # fig, ax = plt.subplots(1)
    # plot_sensor(-5, -20, ax)


    # ax3d = plt.figure().add_subplot(projection="3d")
    # ax3d.plot_surface(np.array([-33.26, -5.1]), np.array([-33.51, 5.1]), np.array([-272.0, -252]))
    # Create axis

    # # # x, y, z = np.indices((35, 40, 200))
    # # print(x,y,z)
    # r_x = [-5.1, 37.41]
    # r_y = [-33.26, -5.1]
    # r_z = [-0.25, 0.25]
    # y = np.arange(-35, 35)
    # x = np.arange(-40, 40)
    # # z = np.arange(np.min(z_left_sensors), np.max(z_right_sensors))
    # z = np.arange(-200, 200)

    # print(combinations(np.array(list(product(r_x, r_y, r_z))), 2))

    # for s, e in combinations(np.array(list(product(r_x, r_y, r_z))), 2):
    #     if np.sum(np.abs(s-e)) == r_x[1]-r_x[0]:
    #         ax3d.plot3D(*zip(s, e), color="green")
    #     elif np.sum(np.abs(s-e)) == r_y[1]-r_y[0]:
    #         ax3d.plot3D(*zip(s, e), color="green")
    #     elif np.sum(np.abs(s-e)) == r_z[1]-r_z[0]:
    #         ax3d.plot3D(*zip(s, e), color="green")

    # plt.show()
