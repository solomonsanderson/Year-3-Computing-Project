'''creating a class for the velo sensor'''


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
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


        for z in self.left_sensor_z:
            l_sensor = left_sensor(z)
            self.l_sensors.append(l_sensor)

        for z in self.right_sensor_z:
            r_sensor = right_sensor(z)
            self.r_sensors.append(r_sensor)


    def hits(self, particle):
        ''''''
        hits = 0

        particle.calibrate(self.zs)

        # particle.calibrate(self.left_sensor_z)  # may be able to do this above
        # particle.calibrate(self.right_sensor_z)

        hits = []

        for i in range(len(self.left_sensor_z)):
            l_x, l_y, l_z = particle.position(self.left_sensor_z[i])
            l_hits = self.l_sensors[i].in_sensor(l_x, l_y)
            print(f"left hits {l_hits}")
            if (l_hits[0].size > 0 ) & (l_hits[1].size > 0):
                hits.append((l_z, l_x, l_y))
            

        for i in range(len(self.right_sensor_z)):
            r_x, r_y, r_z = particle.position(self.right_sensor_z[i])
            r_hits = self.r_sensors[i].in_sensor(r_x, r_y)

            print(f"r hits {r_hits}")
            # print(r_hits[0].size)
            # print(r_hits[1].size)
            # print(r_hits[0].size > 0)
            if (r_hits[0].size > 0):
                hits.append((r_z, r_x, r_y))

            # print(f"right hits {r_hits}")


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
    
    hits = velo(z_left_sensors, z_right_sensors).hits(part)
    print(hits)

    z_all = np.concatenate([z_left_sensors, z_right_sensors])
    # l_x, l_y, l_z = part.position(z_left_sensors)  # may be able to do this above
    # r_x, r_y, r_z = part.position(z_right_sensors)
    x, y, z = part.position(z_all)
    
    ax3d = plt.figure().add_subplot(projection="3d")

    # Right Sensors
    
    positions = [(0, -5.1, -33.26), (0, 5.1, -5.1)]
    positions_arr = []
    
    for i in z_right_sensors:
        # print(f"right {i}")
        positions_arr.append((i, -5.1, -33.26))
        positions_arr.append((i, 5.1, -5.1))

        
    
    sizes = [(0.2, 42.51, 28.41), (0.2, 28.16, 38.61)]
    sizes_arr = sizes * len(z_right_sensors)

    ones = np.ones(len(z_right_sensors))
    colors_arr = [("crimson")] * len(z_right_sensors) * 2

    # print(colors_arr)
    # print(sizes_arr)
    # print(positions_arr)

    pc = plotCubeAt2(positions_arr, sizes_arr, colors=colors_arr, edgecolor="k")
    ax3d.add_collection3d(pc)    


    # Left Sensors
    positions = [(0, -33.06, -33.26),(0, -37.41, 5.1)]  # z, x, y
    # sizes = [(0.2, 28.16, 38.61), (0.2, 42.51, 28.41)]  # z, x, y
    # colors = [("blue"), ("blue")]

    positions_arr = []
    for i in z_left_sensors:
        # print(f"left {i}")
        positions_arr.append((i, -33.06, -33.26))
        positions_arr.append((i, -37.41, 5.1))

    
    sizes = [(0.2, 28.16, 38.61), (0.2, 42.51, 28.41)]
    sizes_arr = sizes * len(z_left_sensors)

    ones = np.ones(len(z_left_sensors))
    colors_arr = [("blue")] * len(z_right_sensors) * 2
    
    # print(len(positions_arr))
    # print(len(sizes_arr))    
    # print(len(colors_arr))


    pc = plotCubeAt2(positions_arr, sizes_arr, colors=colors_arr, edgecolor="k")
    ax3d.add_collection3d(pc)    

    # print(z_right_sensors)
    # print(part.z_arr)
    ax3d.plot(part.z_arr, part.x_arr,  part.y_arr,  marker="o", color="green")
    ax3d.scatter(0,0,0, color="y", s=15)
    # print(hits)
    for hit in hits:
        # print(hit)
        ax3d.scatter(hit[0], hit[1], hit[2], color="hotpink")

    # ax3d.scatter(0, -5.1, -5.1)
    ax3d.set_xlabel("z")
    ax3d.set_ylabel("x")
    ax3d.set_zlabel("y")
    ax3d.set_zlim(-50, 50)
    ax3d.set_ylim(-50,50)
    ax3d.set_xlim(-300, 800)

    legend_elements = [Line2D([0], [0], color="green", lw=4, label="Particle Path"),
                Patch(facecolor='crimson', edgecolor='black',
                         label='Right Sensors'),
                Patch(facecolor='blue', edgecolor='black',
                         label='Left Sensors')]

    ax3d.legend(handles = legend_elements)
    plt.show()