'''

'''


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from particle import particle
from velo import velo, fit_3d, particle_generate, uniform_particle_generate, recon_eff, recon_eta
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
    hits = velo_detect.hits(part)
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

    particles = particle_generate(50, 4, z_all)  # Generatess 50 particles with a pseudo rapidity of 4
    reconstruction_eff, reconstruction_hits = recon_eff(particles, velo_detect)
    len(reconstruction_hits)
    print(f"reconstruction efficiency {reconstruction_eff}")
    # print(len(reconstruction_eff))
    par = particles[0]

    # recon_hits = []
    # for par in particles:
    #     hits = velo_detect.hits(par)
    #     ax3d.plot(par.z_arr, par.x_arr,  par.y_arr,  marker=None, alpha=0.5, color="green")
    #     if len(hits) >=3:
    #         recon_hits.append(hits)
    
    # print(np.array(recon_hits[0])[])
    # print(z, x, y)
    # fit_3d(z, x, y)

    # ax3d.scatter(par.z_arr, par.x_arr,  par.y_arr,  marker=None, alpha=0.5, color="green")
    # linepts = fit_3d(particles[0].z_arr, particles[0].x_arr, particles[0].y_arr)
    # ax3d.plot(*linepts.T, color="orange")
    # points = list(zip(par.z_arr, par.x_arr, par.y_arr))
    # line_fit = Line.  (points)

    
    rand_rapidity = np.random.uniform(low=0.1, high=7.04, size = 100 )
    # rapid_fig, rapid_ax = plt.subplots()
    unif_rapidity = np.linspace(0, 10, 100)
    effs = []
    count = 0

    etas, effs = recon_eta([0, 10], z_all, velo_detect)    
    fig_eff, ax_eff = plt.subplots()
    ax_eff.plot(etas, effs)


    p_ts = []
    sigmas = []
    etas = []
    phis = []
    uniform_particles, phi_values, eta_values = uniform_particle_generate((0, 2 * np.pi), (5, 6), 1000, z_all)
    for count, particle_ in enumerate(uniform_particles):
        particle_.set_pmag(10)  # setting the particles total momentum in GeV
        hits = velo_detect.hits(particle_)
        if len(hits) >= 3:
            # ax3d.plot(particle_.z_arr, particle_.x_arr,  particle_.y_arr,  marker=None, alpha=0.5, color="green")
            # print(np.array(hits)[:,0])
            z, x, y = np.array(hits)[:,0].flatten(), np.array(hits)[:,1].flatten(), np.array(hits)[:,2].flatten()
            fit_result  = fit_3d(z, x, y, plot=False)
            # print(f"fit: {fit_result}")
            p_t = particle_.transverse_momentum(*fit_result[0:2], 10)
            sigma_p_t = particle_.p_resolution(*fit_result)
            
            # print(f"p_t = {p_t}, sigma_p_t = {sigma_p_t}")
            p_ts.append(p_t)
            etas.append(eta_values[count])
            phis.append(phi_values[count])
            sigmas.append(sigma_p_t)
            
    # uniform_particles[3].set_pmag(10)
    
    pt_fig, pt_ax = plt.subplots(2)
    # print(len(phis), len(sigma_p_t))
    # pt_ax[0].errorbar(phis, sigma_p_t)
    # pt_ax[1].errorbar(etas, sigma_p_t)

    n_bins = 24
    width = (np.max(etas) - np.min(etas))/n_bins
    bins=[ np.min(etas) + i*width for i in range(n_bins+1) ]
    bin_centres = [ (bins[i] + bins[i+1])/2 for i in range(len(bins)-1)]
    print(f"bin pos{bin_centres}")

    indices  = np.digitize(etas, bins) #< 100 indices 1
    print(indices)

    pt_plot = [ 0 for i in range(n_bins)]
    print(pt_plot)
    events_per_bin = [ 0 for i in range(len(bins) -1 )]
    for i, index in enumerate(indices):
        # print(i, index)
        pt_plot[index - 2] += p_ts[i]
        events_per_bin[index - 2] += 1

    pt_ax[0].bar(bin_centres, [pt/n if n > 0 else 0 for pt, n in zip(pt_plot,events_per_bin)], width=width, align="center")


    print(sorted(p_ts))
    # pt_ax[0].plot(phis, p_ts)
    # pt_ax[0].hist(p_ts)
    pt_ax[0].set_xlabel("$\eta$")
    pt_ax[0].set_ylabel("$P_T$")
    pt_ax[1].plot(etas, p_ts)
    pt_ax[1].set_xlabel("$\eta$")
    pt_ax[1].set_ylabel("$P_T$")
    pt_fig.tight_layout()
    plt.show()