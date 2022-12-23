''' sensor.py

This script allows the user to create left and right sensor arrays for the VELO
detector at LHCb. 

This script requires numpy and the script cube_plot.py.

This script contains the following functions:
    * plot_sensor - plots one of the velo sensors in 2 dimensions.
    * plot_sensor_3d - plots one of the velo sensors in 3 dimensions.

It also contains the following classes:
    * right_sensor - class describing the right hand sensors of the VELO
      detector. 
    * left_sensor - class describing the left hand sensors of the VELO detector.
'''


import matplotlib.pyplot as plt
import numpy as np
from cube_plot import plotCubeAt2



# We have 26 sensors for each side. 
z_left_sensors = np.array([-277, -252, -227, -202, -132, -62, -37, -12, 13, 38, 63, 88, 
113, 138, 163, 188, 213, 238, 263, 325, 402, 497, 616, 661, 706, 751])

z_left_sensors_bounds = np.array(list(zip(z_left_sensors, z_left_sensors + 0.2)))

z_right_sensors = np.array([-289, -264, -239, -214, -144, -74, -49, -24, 1, 26, 51, 76, 
101, 126, 151, 176, 201, 226, 251, 313, 390, 485, 604, 649, 694, 739])


def plot_sensor(x, y, axis):
    '''Function for making a 2 dimensional plot of the sensors and any hits on
    them.
    
    Args:
        x - float/int, position of particle hits to be plotted on the sensor.
        y - float/int, position of particle hots to be plotted on the sensor.
        axis - matplotlib ax object, the axis for the sensor and the particle
        hit to be plotted on.
    
    Returns:
        None.
    '''

    
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


def plot_sensor_3d(axis):
    '''Function for plotting the sensors in 3 dimensions. Plots both the left
    and right sensors.
    
    Args:
        axis - matplotlib axis object with 3d projection, the axis for the
        sensor and particle track to be plotted on.

    Returns: 
        None.
    '''


    # Right Sensors
    positions_arr = []
    # making an array of positions for right sensors.
    for i in z_right_sensors:
        positions_arr.append((i, -5.1, -33.26))
        positions_arr.append((i, 5.1, -5.1))

    sizes = [(0.2, 42.51, 28.41), (0.2, 28.16, 38.61)]  # settings sizes of sensors
    sizes_arr = sizes * len(z_right_sensors)
    colors_arr = [("crimson")] * len(z_right_sensors) * 2  # making an array of colors for each rectangle

    pc = plotCubeAt2(positions_arr, sizes_arr, colors=colors_arr, edgecolor="k")
    axis.add_collection3d(pc)    


    # Left Sensors
    positions_arr = []

    # making an array of positions for left sensors.
    for i in z_left_sensors:
        positions_arr.append((i, -33.06, -33.26))
        positions_arr.append((i, -37.41, 5.1))
  
    sizes = [(0.2, 28.16, 38.61), (0.2, 42.51, 28.41)]
    sizes_arr = sizes * len(z_left_sensors)
    colors_arr = [("blue")] * len(z_right_sensors) * 2  # making an array of colors for each rectangle

    pc = plotCubeAt2(positions_arr, sizes_arr, colors=colors_arr, edgecolor="k")  # creating cuboids
    axis.add_collection3d(pc)  # adding objects to the axis.
    axis.scatter(0,0,0, color="black", s=15)  # sets point at the origin
    # Formatting plot
    axis.set_xlabel("z")
    axis.set_ylabel("x")
    axis.set_zlabel("y")
    axis.set_zlim(-50, 50)
    axis.set_ylim(-50,50)
    axis.set_xlim(-300, 800)


class right_sensor():
    ''' A class defining a right hand sensor in the vertex locator at LHCb.
    
    Attributes:
        z_pos - array/list, the z posions of the right sensors.'''

 
    def __init__(self, z_pos):
        '''Constructor for right_sensor.
        Args:
            z_pos - array/list, the z posions of the right sensors.
        Returns:
            None.
         '''
        self.z_pos = z_pos


    def in_sensor(self, x, y):
        '''Calculating if the particle lies within the right hand sensor.
        
        Args: 
            x: float - x position of particle.
            y: float - y position of particle.
        
        Returns: 
            (hits_y, hits_x) - tuple, a tuple of arrays
            of the y and x positions of particle hits.
            '''


        index = np.where((((5.1 < x) & (x < 33.26)) & ((-5.1 < y) & (y < 33.51))) | (((-5.1 < x) & (x < 37.41)) & ((-33.51 < y) & (y < -5.1))))  # get the indicies where the particle is inside the sensor.
        hits_y, hits_x = y[index], x[index] # return x and y values when the particle is inside
        return hits_y, hits_x


class left_sensor():
    ''' A class defining a left hand sensor in the vertex locator at LHCb.
    
    Attributes:
        z_pos - array/list, the z posions of the left sensors.'''
        


    def __init__(self, z_pos):
        '''Constructor for right_sensor.
        Args:
            z_pos - array/list, the z posions of the right sensors.
        Returns:
            None.
         '''

        
        self.z_pos = z_pos


    def in_sensor(self, x, y):
        '''Calculating if the particle lies within the left hand sensor.
        
        Args: 
            x: float - x position of particle.
            y: float - y position of particle.
        
        Returns: 
            (hits_y, hits_x) - tuple, a tuple of arrays
            of the y and x positions of particle hits.
            '''
        

        index = np.where((((-33.26 < x) & (x < -5.1)) & ((-33.51 < y) & (y < 5.1))) | (((-37.41 < x) & (x < 5.1)) & ((5.1 < y) & (y < 33.51))))  # get the indicies where the particle is inside the sensor.
        hits_y, hits_x = y[index], x[index]  # return x and y values when the particle is inside
        
        return hits_y, hits_x
        

if __name__ == "__main__":
    pass