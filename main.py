'''main.py
The main file of my project where all of the other files are called from and
used. 

This file outputs...
'''


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

from particle import particle
from velo import velo, fit_3d
from sensor import plot_sensor_3d



if __name__ == "__main__":
    #  defining arrays of z positions for both left and right sensors
    z_right_sensors = np.array([-289, -264, -239, -214, -144, -74, -49, -24, 1, 26, 51, 76, 
    101, 126, 151, 176, 201, 226, 251, 313, 390, 485, 604, 649, 694, 739])

    z_left_sensors = np.array([-277, -252, -227, -202, -132, -62, -37, -12, 13, 38, 63, 88, 
    113, 138, 163, 188, 213, 238, 263, 325, 402, 497, 616, 661, 706, 751])

    
    velo_detect = velo(z_left_sensors, z_right_sensors)
    z_all = np.concatenate([z_left_sensors, z_right_sensors])

    ax3d = plt.figure().add_subplot(projection="3d")
    plot_sensor_3d(ax3d)

    # Formatting our plot including legend
    ax3d.scatter(0,0,0, color="y", s=15)  # places a point at the origin

    ax3d.set_xlabel("z")
    ax3d.set_ylabel("x")
    ax3d.set_zlabel("y")
    ax3d.set_zlim(-50, 50)
    ax3d.set_ylim(-50,50)
    ax3d.set_xlim(-400, 800)


    #  Formatting legend
    legend_elements = [Line2D([0], [0], color="green", lw=4, label="Particle Path"),
                Patch(facecolor='crimson', edgecolor='black',
                           label='Right Sensors'),
                Patch(facecolor='blue', edgecolor='black',
                         label='Left Sensors')]

    ax3d.legend(handles = legend_elements)

    particles = velo_detect.uniform_particle_generate((0, 2 * np.pi), (4,4), 50)[0]  # Generatess 50 particles with a pseudo rapidity of 4
    reconstruction_eff, reconstruction_hits = velo_detect.recon_eff()
    print(f"reconstruction efficiency {reconstruction_eff}")



    zero_etas, zero_effs = velo_detect.recon_eta([-12, 12], ([0, 0], [0, 0], [0, 0]))
    etas, effs = velo_detect.recon_eta([-12, 12], ([-5, 5], [-5, 5], [0, 0]))

    
    fig_eff, ax_eff = plt.subplots()
    ax_eff.plot(zero_etas, zero_effs, color="blue", alpha = 0.5)
    ax_eff.plot(etas, effs, color="red", alpha = 0.5)
    
    ax_eff.set_xlabel("Pseudorapdidity, $\eta$")
    ax_eff.set_ylabel("Reconstruction Efficiency")
    ax_eff.set_title("Plot of Track Reconstruction Efficiency vs. Pseudorapidity")

    true_p_ts = []
    p_ts = []
    sigmas = []
    etas = []
    phis = []
    uniform_particles, phi_values, eta_values = velo_detect.uniform_particle_generate((0, 2  * np.pi), (-6, 6), 1000, ([-5, 5], [-5, 5], [0, 0]))
    for count, particle_ in enumerate(uniform_particles):
        particle_.set_pmag(10)  # setting the particles total momentum in GeV
        hits = velo_detect.hits(particle_, hit_resolution=5)
        if len(hits) >= 3:
            # ax3d.plot(particle_.z_arr, particle_.x_arr,  particle_.y_arr,  marker=None, alpha=0.5, color="green")
            z, x, y = np.array(hits)[:,0].flatten(), np.array(hits)[:,1].flatten(), np.array(hits)[:,2].flatten()
            fit_result  = fit_3d(z, x, y, plot=False)
            p_t = particle_.transverse_momentum(*fit_result[0:2], 10)
            sigma_p_t = particle_.p_resolution(*fit_result[0:4])
            true_p_t = (particle_.p_arr[0] ** 2 +  particle_.p_arr[1] ** 2 ) ** 0.5  # calculating the momentum from the initial momentum values.
            true_p_ts.append(true_p_t)
            # print(f"x{particle_.p_x}, y{particle_.p_y}")
            # print(f"p_t = {p_t}, sigma_p_t = {sigma_p_t}, true p_t = {true_p_t}")
            p_ts.append(p_t)
            etas.append(eta_values[count])
            phis.append(phi_values[count])
            sigmas.append(sigma_p_t)
            
    
    sigmas = np.array(sigmas)
    p_ts = np.array(p_ts)
    hist_fig, hist_ax = plt.subplots(1)
    # print(sigmas)
    # print(np.min(sigmas/p_ts), np.max(sigmas/p_ts))
    print(f"momentum mean = {np.mean(p_ts)} +/- {np.mean(sigmas)}")
    hist_ax.hist(sigmas/p_ts,bins = 100 , range = (-0.05, 0.05))
    pt_fig, pt_ax = plt.subplots(2)

    eta_lo = -6
    eta_hi = 6
    n_bins = 96

    width = (eta_hi - eta_lo)/n_bins
    bins = [eta_lo + i*width for i in range(n_bins + 1)]
    bin_centres = [(bins[i] + bins[i+1])/2 for i in range(len(bins) -1 )]

    indices  = np.digitize(etas, bins) #< 100 indices 1
    pt_plot = np.zeros(n_bins)
    events_per_bin = np.zeros(n_bins, dtype=float)
    true_pt_plot = np.zeros(n_bins)
    true_events_per_bin = np.zeros(n_bins, dtype=float)


    for i, index in enumerate(indices):
        pt_plot[index - 2] += p_ts[i]
        events_per_bin[index - 2] += 1
        true_pt_plot[index - 2] += true_p_ts[i]
        true_events_per_bin[index - 2] += 1
    
    
    pt_ax[0].bar(bin_centres, pt_plot / events_per_bin, width=width, align="center", alpha=0.5, label = "fitted")
    pt_ax[0].bar(bin_centres, true_pt_plot / true_events_per_bin, width=width, align="center", alpha=0.5, label="true")

    pt_ax[0].set_ylabel("$P_T$")
    pt_ax[0].set_xlabel("$\eta$")
    pt_ax[0].legend()
    pt_fig.tight_layout()
    pt_ax[0].set_title("Plot of $P_T$ against $\eta$. With 10000 Points")

    
    
    
    phi_lo = 0
    phi_hi = 2 * np.pi
    n_bins = 96

    width = (phi_hi - phi_lo)/n_bins
    bins = [phi_lo + i*width for i in range(n_bins + 1)]
    bin_centres = [(bins[i] + bins[i+1])/2 for i in range(len(bins) -1 )]
    # print(f"bin cents {bin_centres}")
    # print(f"bin pos{bin_centres}")

    indices  = np.digitize(phis, bins) #< 100 indices 1
    pt_plot = np.zeros(n_bins)
    events_per_bin = np.zeros(n_bins, dtype=float)


    for i, index in enumerate(indices):
        pt_plot[index - 2] += p_ts[i]
        events_per_bin[index - 2] += 1


    pt_ax[1].bar(bin_centres, pt_plot / events_per_bin, width=width, align="center")
    

    # print(sorted(p_ts))
    # pt_ax[0].plot(phis, p_ts)
    # pt_ax[0].hist(p_ts)
    pt_ax[1].set_ylabel("$P_T$")
    pt_ax[1].set_xlabel("$\phi$")
    # pt_ax[1].plot(etas, p_ts)
    # pt_ax[1].set_xlabel("$\eta$")
    # pt_ax[1].set_ylabel("$P_T$")
    pt_fig.tight_layout()
    pt_ax[1].set_title("Plot of $P_T$ against $\phi$. With 10000 Points")

    # investigating the effect of different starting positions on the reconstruction efficiency.
    print("range of starts")
    # velo_detector = velo(z_left_sensors, z_right_sensors)
    # particles = velo_detector.uniform_particle_generate((0, 2 * np.pi), (4,4), 50, ([0,0],[0,0],[0, 0]))[0]  # Generatess 50 particles with a pseudo rapidity of 4
    reconstruction_eff, reconstruction_hits = velo_detect.recon_eff()
    print(f"range reconstruction efficiency {reconstruction_eff}")

    # print(len(velo_detect.particles))
    velo_detect.uniform_particle_generate((0, 2  * np.pi), (-6, 6), 100, ([-10, 10], [-10, 10], [0, 0]))
    ips = []
    sigma_ips = []
    for par in velo_detect.particles:
        ax3d.plot(par.z_arr, par.x_arr,  par.y_arr,  marker=None, alpha=0.5, color="green")
        ips.append(par.impact_parameter(*fit_result[0:2], *fit_result[4:6]))
        sigma_ips.append(par.ip_resolution(*fit_result))
    # print(ips)
    # print(np.mean(ips), np.mean(sigma_ips))

    # fig_zvar, ax_zvar = plt.subplots()
    # ax_zvar.set_xlabel("Pseudorapdidity, $\eta$")
    # ax_zvar.set_ylabel("Reconstruction Efficiency")
    # ax_zvar.set_title("Plot of Track Reconstruction Efficiency vs. Pseudorapidity")
    plt.show()