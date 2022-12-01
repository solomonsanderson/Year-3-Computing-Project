'''creating a class for the velo sensor'''


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from sensor import left_sensor, right_sensor
from particle import particle
from sensor import plot_sensor, plot_sensor_3d
from itertools import product, combinations
import random
from scipy.optimize import curve_fit
# from skspatial.objects import Line, Points
# from skspatial.plotting import plot_3d



def pr_theta(pseudo_rapidity):
    ''''''
    theta = 2 * np.arctan(np.exp(-pseudo_rapidity))
    return theta


def particle_generate(number, pseudo_rapidity, z):
    theta = pr_theta(pseudo_rapidity)
    ys = z * np.tan(theta)

    
    particles = []


    # xys = list(zip(x,y))
    while len(particles) < number:
        phi = np.random.uniform(0, 2 * np.pi)
        xs = z * ((np.tan(theta))/(np.tan(phi)))
        # print(xy, z)
        # print(*xy, z[count])
        part_obj = particle((xs[0] ,ys[0], z[0]), (0,0,0))
        part_obj.calibrate(z)
        particles.append(part_obj)
    
    return particles


def recon_eff(particles):
    reconned = 0
    for par in particles:
        # print(par)
        hits, number = velo_detect.hits(par)
        if number >= 3:
            reconned += 1

        # print(par.z_ar/r, par.x_arr,  par.y_arr)
        ax3d.plot(par.z_arr, par.x_arr,  par.y_arr,  marker=None)

    print(reconned)
    print(len(particles))
    efficiency = reconned / len(particles)
    return efficiency


def fit_3d(x, y, z):
    data = np.array(list(zip(x,y,z)))
    mean = data.mean(axis=0)
    # print(mean)
    uu, dd, vv = np.linalg.svd(data - mean)
    first_pc = vv[0]
    linepts = vv[0] * np.mgrid[-400:800:1000j][:, np.newaxis]
    print(linepts)
    linepts += mean
    return linepts



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


    def hits(self, particle, hit_efficiency = 98, hit_resolution = 0.012):
        ''''''
        hits = 0

        particle.calibrate(self.zs)
        
        # particle.calibrate(self.left_sensor_z)  # may be able to do this above
        # particle.calibrate(self.right_sensor_z)

        hits = []

        for i in range(len(self.left_sensor_z)):
            l_x, l_y, l_z = particle.position(self.left_sensor_z[i])
            l_x_smear = np.random.normal(l_x, hit_resolution)
            l_y_smear = np.random.normal(l_y, hit_resolution)

            # If we dont use the smeared values for l_hits there is no risk of
            # us losing or gaining values through smearing them out of the sensor
            l_hits = self.l_sensors[i].in_sensor(l_x, l_y)


            if (l_hits[0].size > 0 ) & (l_hits[1].size > 0) & (random.choices([True, False], weights = (hit_efficiency, 100 - hit_efficiency))[0]):
                hits.append((l_z, l_x_smear, l_y_smear))
            

        for i in range(len(self.right_sensor_z)):
            r_x, r_y, r_z = particle.position(self.right_sensor_z[i])
            r_x_smear = np.random.normal(r_x, hit_resolution)
            r_y_smear = np.random.normal(r_y, hit_resolution)
            r_hits = self.r_sensors[i].in_sensor(r_x, r_y)


            if (r_hits[0].size > 0) & (random.choices([True, False], weights = (hit_efficiency, 100 - hit_efficiency))[0]):
                hits.append((r_z, r_x_smear, r_y_smear))


        number = len(hits)

        return hits, number


if __name__ == "__main__":
    z_right_sensors = np.array([-289, -264, -239, -214, -144, -74, -49, -24, 1, 26, 51, 76, 
    101, 126, 151, 176, 201, 226, 251, 313, 390, 485, 604, 649, 694, 739])

    z_left_sensors = np.array([-277, -252, -227, -202, -132, -62, -37, -12, 13, 38, 63, 88, 
    113, 138, 163, 188, 213, 238, 263, 325, 402, 497, 616, 661, 706, 751])

    # generating a random phi value:
    phi = np.random.uniform(0, 2 * np.pi)

    part = particle((np.random.uniform(-1,1), np.random.uniform(-1, 1), np.random.uniform(-1, 1)), (0, 0, 0))

    
    velo_detect = velo(z_left_sensors, z_right_sensors)
    hits, number = velo_detect.hits(part)
    # print(hits)

    z_all = np.concatenate([z_left_sensors, z_right_sensors])

    x, y, z = part.position(z_all)
    
    ax3d = plt.figure().add_subplot(projection="3d")
    plot_sensor_3d(ax3d)


    # print(z_right_sensors)
    # print(part.z_arr)
    # ax3d.plot(part.z_arr, part.x_arr,  part.y_arr,  marker="o", color="green")
    ax3d.scatter(0,0,0, color="y", s=15)
    # print(hits)
    # for hit in hits:
    #     # print(hit)
    #     ax3d.scatter(hit[0], hit[1], hit[2], color="hotpink")

    # ax3d.scatter(0, -5.1, -5.1)
    ax3d.set_xlabel("z")
    ax3d.set_ylabel("x")
    ax3d.set_zlabel("y")
    ax3d.set_zlim(-50, 50)
    ax3d.set_ylim(-50,50)
    ax3d.set_xlim(-400, 800)


    legend_elements = [Line2D([0], [0], color="green", lw=4, label="Particle Path"),
                Patch(facecolor='crimson', edgecolor='black',
                           label='Right Sensors'),
                Patch(facecolor='blue', edgecolor='black',
                         label='Left Sensors')]

    ax3d.legend(handles = legend_elements)

    particles = particle_generate(1, 4, z_all)
    print(f" reconstruction efficiency {recon_eff(particles)}")
    par = particles[0]

    # for par in particles:
    #     hits, number = velo_detect.hits(par)
    #     ax3d.plot(par.z_arr, par.x_arr,  par.y_arr,  marker=None, alpha=0.5, color="green")
    ax3d.scatter(par.z_arr, par.x_arr,  par.y_arr,  marker=None, alpha=0.5, color="green")
    linepts = fit_3d(particles[0].z_arr, particles[0].x_arr, particles[0].y_arr)
    ax3d.plot(*linepts.T, color="orange")
    # points = list(zip(par.z_arr, par.x_arr, par.y_arr))
    # line_fit = Line.best_fit(points)

    plt.show()
