'''

'''


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from particle import particle
from velo import velo, fit_3d, particle_generate, uniform_particle_generate, recon_eff
from sensor import plot_sensor_3d





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

    particles = particle_generate(50, 4, z_all)
    reconstruction_eff, reconstruction_hits = recon_eff(particles, velo_detect)
    len(reconstruction_hits)
    print(f"reconstruction efficiency {reconstruction_eff}")
    # print(len(reconstruction_eff))
    par = particles[0]

    for par in particles:
        hits, number = velo_detect.hits(par)
        ax3d.plot(par.z_arr, par.x_arr,  par.y_arr,  marker=None, alpha=0.5, color="green")
    # ax3d.scatter(par.z_arr, par.x_arr,  par.y_arr,  marker=None, alpha=0.5, color="green")
    linepts = fit_3d(particles[0].z_arr, particles[0].x_arr, particles[0].y_arr)
    # ax3d.plot(*linepts.T, color="orange")
    # points = list(zip(par.z_arr, par.x_arr, par.y_arr))
    # line_fit = Line.best_fit(points)

    
    rand_rapidity = np.random.uniform(low=0.1, high=7.04, size = 100 )
    rapid_fig, rapid_ax = plt.subplots()
    unif_rapidity = np.linspace(0, 10, 100)
    effs = []
    count = 0
    for num in unif_rapidity:
        count +=1
        print(count)
        particles = particle_generate(10, num, z_all)
        recon_efficiency = recon_eff(particles, velo_detect)[0]
        effs.append(recon_efficiency)
    
    fig, ax = plt.subplots(1)
    ax.plot(unif_rapidity, effs)
    ax.set_xlabel("Rapidity $\eta$")
    ax.set_ylabel("Reconstruction efficiency")


    # print(linepts[0])
    # print(linepts[1])
    # slope = np.sqrt((linepts[1][0] - linepts[0][0]) ** 2 + (linepts[1][1] - linepts[0][1]) ** 2 + (linepts[1][2] - linepts[0][2]) ** 2)
    # print(slope)
    # print(list(zip(rand_rapidity, effs)))
    

    # uniform sampling
    uniform_particles = uniform_particle_generate((0, 2 * np.pi), (3.5, 4.5), 10, z_all)
    square = uniform_particles[0].x_arr**2 + uniform_particles[0].y_arr**2 + uniform_particles[0].z_arr**2
    root = square ** 0.5
    print(root)
    plt.show()