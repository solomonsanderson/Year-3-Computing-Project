'''main.py

!RUN THE CODE FROM THIS FILE!

The main file of my project where all of the other files are called from and
used. 

This file outputs the following:
    - A plot of reconsturction efficiency for particles generated at the origin
    vs those generated over random x y and z values.
    - Values for the transverse momentum of our particle and its resolution.
    - The impact parameter and its resolution.
    - A 3d plot of the particles and the sensor.
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

    velo_detect = velo(z_left_sensors, z_right_sensors) # creating a velo object
    z_all = np.concatenate([z_left_sensors, z_right_sensors])

    np.seterr(divide="ignore")

    # Plot of particles with varied etas
    ax3d_eta = plt.figure().add_subplot(projection="3d")
    plot_sensor_3d(ax3d_eta)

    lo_eta_particles, lo_eta_values, lo_eta_values = velo_detect.uniform_particle_generate((0, 2  * np.pi), (1, 2), 50 , ([-15e-3, 15e-3], [-15e-3, 15e-3], [-63, 63]))
    for par in lo_eta_particles:
        ax3d_eta.plot(par.z_arr, par.x_arr,  par.y_arr,  marker=None, alpha=0.5, color="green")

    hi_eta_particles, hi_eta_values, hi_eta_values = velo_detect.uniform_particle_generate((0, 2  * np.pi), (5, 6), 50 , ([-15e-3, 15e-3], [-15e-3, 15e-3], [-63, 63]))
    for par in hi_eta_particles:
        ax3d_eta.plot(par.z_arr, par.x_arr,  par.y_arr,  marker=None, alpha=0.5, color="orange")


    ax3d_eta.set_title("Plot of particle tracks at high and low $\eta$.")
    ax3d_eta.set_xlabel("z"), ax3d_eta.set_ylabel("x"), ax3d_eta.set_zlabel("y")
    ax3d_eta.set_zlim(-50, 50), ax3d_eta.set_ylim(-50,50), ax3d_eta.set_xlim(-400, 800)

    eta_legend_elements = [Line2D([0], [0], color="green", lw=4, label="$1\\leq \\eta \\leq 2$"), 
                Line2D([0], [0], color="orange", lw=4, label="$5\\leq \\eta \\leq 6$"),
                Patch(facecolor='crimson', edgecolor='black',
                           label='Right Sensors'),
                Patch(facecolor='blue', edgecolor='black',
                         label='Left Sensors')]

    ax3d_eta.legend(handles = eta_legend_elements)

    # Plot of particles with varied phis
    ax3d_phi = plt.figure().add_subplot(projection="3d")
    plot_sensor_3d(ax3d_phi)
    lo_phi_particles, lo_phi_values, lo_phi_values = velo_detect.uniform_particle_generate((0, 0.5  * np.pi), (-6, 6), 50 , ([-15e-3, 15e-3], [-15e-3, 15e-3], [-63, 63]))
    for par in lo_phi_particles:
        ax3d_phi.plot(par.z_arr, par.x_arr,  par.y_arr,  marker=None, alpha=0.5, color="green")

    hi_phi_particles, hi_phi_values, hi_phi_values = velo_detect.uniform_particle_generate((1.5 * np.pi, 2  * np.pi), (-6, 6), 50 , ([-15e-3, 15e-3], [-15e-3, 15e-3], [-63, 63]))
    for par in hi_phi_particles:
        ax3d_phi.plot(par.z_arr, par.x_arr,  par.y_arr,  marker=None, alpha=0.5, color="orange")


    ax3d_phi.set_title("Plot of particle tracks at varied $\phi$ values.")
    ax3d_phi.set_xlabel("z"), ax3d_phi.set_ylabel("x"), ax3d_phi.set_zlabel("y")
    ax3d_phi.set_zlim(-50, 50), ax3d_phi.set_ylim(-50,50), ax3d_phi.set_xlim(-400, 800)

    phi_legend_elements = [Line2D([0], [0], color="green", lw=4, label="$0\\leq \\phi \\leq \\frac{1}{2}\\pi$"), 
                Line2D([0], [0], color="orange", lw=4, label="$\\frac{3}{2} \\pi \\leq \\phi \\leq 2 \\pi$"),
                Patch(facecolor='crimson', edgecolor='black',
                           label='Right Sensors'),
                Patch(facecolor='blue', edgecolor='black',
                         label='Left Sensors')]

    ax3d_phi.legend(handles = phi_legend_elements)

    # reconstruction efficiency as a function of pseudorapidity
    zero_etas, zero_effs = velo_detect.recon_eta([-7.5, 7.5], ([0, 0], [0, 0], [0, 0]))  # origin
    etas, effs = velo_detect.recon_eta([-7.5, 7.5], ([-15e-3, 15e-3], [-15e-3, 15e-3], [-60, 60]))  # region

    # Plotting reconstruction efficiencies
    fig_eff, ax_eff = plt.subplots()
    ax_eff.plot(zero_etas, zero_effs
    , color="blue", alpha = 0.5, label="Origin")
    ax_eff.plot(etas, effs, color="red", alpha = 0.5, label="Random")
    # Formatting plot
    ax_eff.legend()
    ax_eff.set_xlabel("Pseudorapdidity, $\eta$")
    ax_eff.set_ylabel("Reconstruction Efficiency")
    ax_eff.set_title("Plot of Track Reconstruction Efficiency vs. Pseudorapidity")

    print("Reconstruction Efficiencies Calculated")

    true_p_ts = []
    p_ts = []
    sigmas = []
    etas = []
    phis = []
    uniform_particles, phi_values, eta_values = velo_detect.uniform_particle_generate((0, 2  * np.pi), (-7.5, 7.5), 10000 , ([-15e-3, 15e-3], [-15e-3, 15e-3], [-63, 63]))
    for count, particle_ in enumerate(uniform_particles):
        particle_.set_pmag(10)  # setting the particles total momentum in GeV
        hits = velo_detect.hits(particle_)  # getting coordinates of particle hits
        if len(hits) >= 3:
            z, x, y = np.array(hits)[:,0].flatten(), np.array(hits)[:,1].flatten(), np.array(hits)[:,2].flatten()
            fit_result  = fit_3d(z, x, y, plot=False)  # fitting a straight line to the hits

            p_t = particle_.transverse_momentum(*fit_result[0:2], 10)  # calculating the transverse momentum
            sigma_p_t = particle_.p_resolution(*fit_result[0:4])  # calculating the resolution of p_t
            true_p_t = (particle_.p_arr[0] ** 2 +  particle_.p_arr[1] ** 2 ) ** 0.5  # calculating the momentum from the initial momentum values.
            true_p_ts.append(true_p_t)

            p_ts.append(p_t)
            etas.append(eta_values[count])
            phis.append(phi_values[count])
            sigmas.append(sigma_p_t)
            
    
    sigmas = np.array(sigmas)
    p_ts = np.array(p_ts)
    hist_fig, hist_ax = plt.subplots(1)

    print(f"momentum mean = {np.mean(p_ts)} +/- {np.mean(sigmas)}")

    # histogram of the error on p_t
    hist_ax.hist(sigmas/p_ts, bins = 100 , range=(-1e-5, 1e-5))
    hist_ax.set_xlabel("Error on $p_T$")
    hist_ax.set_ylabel("Count")
    hist_ax.set_title("Histogram of the error on $p_T$.")

    pt_fig, pt_ax = plt.subplots(2)

    # bar chart of p_t vs eta
    eta_lo = np.min(etas)
    eta_hi = np.max(etas)
    n_bins = 96

    width = (eta_hi - eta_lo)/n_bins  # calculating bin width
    bins = [eta_lo + i*width for i in range(n_bins + 1)]
    bin_centres = [(bins[i] + bins[i+1])/2 for i in range(len(bins) -1 )]  # calculating positions of bin centres

    indices  = np.digitize(etas, bins)  # returns the indices of the bins in which every input belongs
    pt_plot = np.zeros(n_bins) 
    events_per_bin = np.zeros(n_bins, dtype=float)
    true_pt_plot = np.zeros(n_bins)
    true_events_per_bin = np.zeros(n_bins, dtype=float)

    #  Counting the events in each bin
    for i, index in enumerate(indices):
        pt_plot[index - 2] += p_ts[i]
        events_per_bin[index - 2] += 1
        true_pt_plot[index - 2] += true_p_ts[i]
        true_events_per_bin[index - 2] += 1
    
    
    pt_ax[0].bar(bin_centres, pt_plot / events_per_bin, width=width, align="center", alpha=0.5, label = "fitted", color="red")  # fitted p_t
    pt_ax[0].bar(bin_centres, true_pt_plot / true_events_per_bin, width=width, align="center", alpha=0.5, label="true", color="blue")  # generated p_t

    # formatting plot
    pt_ax[0].set_ylabel("$\sigma_{p_T}/p_T$")
    pt_ax[0].set_xlabel("$\eta$")
    pt_ax[0].legend()
    pt_fig.tight_layout()
    pt_ax[0].set_title("Difference between Fitted and True $p_T$")
    
     # plotting the difference between fitted p_t and true p_t
    pt_ax[1].bar(bin_centres,((pt_plot/events_per_bin) - (true_pt_plot/true_events_per_bin)), width=width ,color="green")
    pt_ax[1].set_ylabel("$\Delta p_T$")
    pt_ax[1].set_xlabel("$\eta$")
    pt_fig.tight_layout()
    pt_ax[1].set_title("Difference between Fitted and True $p_T$.")

    # Plotting 
    pt_phi_fig, pt_phi_ax = plt.subplots()
    phi_lo = 0
    phi_hi = 2 * np.pi
    n_bins = 96

    width = (phi_hi - phi_lo)/n_bins  # calculating bin width
    bins = [phi_lo + i*width for i in range(n_bins + 1)]  
    bin_centres = [(bins[i] + bins[i+1])/2 for i in range(len(bins) -1 )]  # calculating position of bins


    indices  = np.digitize(phis, bins)
    pt_plot = np.zeros(n_bins)
    events_per_bin = np.zeros(n_bins, dtype=float)


    for i, index in enumerate(indices):
        pt_plot[index - 2] += p_ts[i]
        events_per_bin[index - 2] += 1

    
    pt_phi_ax.bar(bin_centres, pt_plot / events_per_bin, width=width, align="center")

    pt_phi_ax.set_ylabel("$P_T$")
    pt_phi_ax.set_xlabel("$\phi$")
    pt_phi_fig.tight_layout()
    pt_phi_ax.set_title("Plot of $P_T$ against $\phi$. With 10000 Points")
   


    # investigating the effect of different starting positions on the reconstruction efficiency.
    reconstruction_eff, reconstruction_hits = velo_detect.recon_eff()
    print(f"range reconstruction efficiency {reconstruction_eff}")

    # print(len(velo_detect.particles))
    velo_detect.uniform_particle_generate((0, 2  * np.pi), (-7.5, 7.5), 10000 , ([0, 0], [0, 0], [0, 0]))
    ips = []
    sigma_ips = []
    
    for par in velo_detect.particles[2:]:
        hits = velo_detect.hits(par)
        if len(hits)>=3:
            z, x, y = np.array(hits)[:,0].flatten(), np.array(hits)[:,1].flatten(), np.array(hits)[:,2].flatten()
            fit_result  = fit_3d(z, x, y, plot=False)
            ips.append(par.impact_parameter(*fit_result[0:2], *fit_result[4:6]))
            sigma_ips.append(par.ip_resolution(*fit_result))


    print(f"IP mean {np.mean(ips)} +/- {np.mean(sigma_ips)}")

    plt.show()