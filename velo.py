'''creating a class for the velo sensor'''

import numpy as np
import matplotlib.pyplot as plt
from sensor import left_sensor, right_sensor
from particle import particle
from sensor import plot_sensor
from itertools import product, combinations
from cube_plot import plotCubeAt2


class velo:
    def __init__(self, left_sensor_z, right_sensor_z):
        self.l_sensors = []
        self.r_sensors = []
        self.left_sensor_z = left_sensor_z
        self.right_sensor_z = right_sensor_z
        self.zs = np.concatenate([left_sensor_z, right_sensor_z])


        for z in range(len(left_sensor_z)):
            l_sensor = left_sensor(z)
            r_sensor = right_sensor(z)
            self.l_sensors.append(l_sensor)
            self.r_sensors.append(r_sensor)


    def hits(self, particle):
        ''''''
        hits = 0

        l_x, l_y = particle.position(self.left_sensor_z)  # may be able to do this above
        r_x, r_y = particle.position(self.right_sensor_z)
        print(l_x)
        fig, ax = plt.subplots(1)
        plot_sensor(l_x, l_y, ax)
        plot_sensor(r_x, r_y, ax)
        plt.show() 

        for i in range(len(self.left_sensor_z)):
            
            if (self.l_sensors[i].in_sensor(l_x, l_y) & self.r_sensors[i].in_sensor(r_x, r_y)) == 0:
                print("error")
            
            elif self.l_sensors[i].in_sensor(l_x, l_y) == 1:
                hits += 1

            elif self.r_sensors[i].in_sensor(r_x, r_y) == 1:
                hits += 1

        return hits


if __name__ == "__main__":
    z_right_sensors = np.array([-289, -264, -239, -214, -144, -74, -49, -24, 1, 26, 51, 76, 
    101, 126, 151, 176, 201, 226, 251, 313, 390, 485, 604, 649, 694, 739])

    z_left_sensors = np.array([-277, -252, -227, -202, -132, -62, -37, -12, 13, 38, 63, 88, 
    113, 138, 163, 188, 213, 238, 263, 325, 402, 497, 616, 661, 706, 751])

    part = particle((np.random.uniform(-1,1), np.random.uniform(-1, 1), np.random.uniform(-1, 1)), (0, 0, 0))
    # part = particle((1,2,3),(0,0,0))

    
    # l_x, l_y = particle.position(self.left_sensor_z)  # may be able to do this above
    # r_x, r_y = particle.position(self.right_sensor_z)
    
    # velo(z_left_sensors, z_right_sensors).hits(part)

    l_x, l_y = part.position(z_left_sensors)  # may be able to do this above
    r_x, r_y = part.position(z_right_sensors)

    ax3d = plt.figure().add_subplot(projection="3d")

    # Right Sensors
    
    positions = [(0, -5.1, -33.26), (0, 5.1, -5.1)]
    positions_arr = []
    
    for i in z_right_sensors:
        positions_arr.append((i, -5.1, -33.26))
        positions_arr.append((i, 5.1, -5.1))

        
    
    sizes = [(0.2, 42.51, 28.41), (0.2, 28.16, 38.61)]
    sizes_arr = sizes * len(z_right_sensors)

    ones = np.ones(len(z_right_sensors))
    colors_arr = [("crimson")] * len(z_right_sensors) * 2

    print(colors_arr)
    print(sizes_arr)
    print(positions_arr)

    pc = plotCubeAt2(positions_arr, sizes_arr, colors=colors_arr, edgecolor="k")
    ax3d.add_collection3d(pc)    


    # Left Sensors
    positions = [(0, -33.06, -33.26),(0, -37.41, 5.1)]  # z, x, y
    # sizes = [(0.2, 28.16, 38.61), (0.2, 42.51, 28.41)]  # z, x, y
    # colors = [("blue"), ("blue")]

    positions_arr = []
    for i in z_right_sensors:
        positions_arr.append((i, -33.06, -33.26))
        positions_arr.append((i, -37.41, 5.1))

    
    sizes = [(0.2, 28.16, 38.61), (0.2, 42.51, 28.41)]
    sizes_arr = sizes * len(z_left_sensors)

    ones = np.ones(len(z_left_sensors))
    colors_arr = [("blue")] * len(z_right_sensors) * 2
    
    print(len(positions_arr))
    print(len(sizes_arr))    
    print(len(colors_arr))


    pc = plotCubeAt2(positions_arr, sizes_arr, colors=colors_arr, edgecolor="k")
    ax3d.add_collection3d(pc)    


    ax3d.plot(z_right_sensors, l_x,  r_y,  marker="o")

    # ax3d.scatter(0, -5.1, -5.1)
    ax3d.set_xlabel("z")
    ax3d.set_ylabel("x")
    ax3d.set_zlabel("y")
    ax3d.set_zlim(-50, 50)
    ax3d.set_ylim(-50,50)
    ax3d.set_xlim(-300, 800)
    plt.show()