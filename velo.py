'''velo.py

Creating a class for the velo sensor using the left_sensor and right_sensor
classes from sensor.py.

This script requires the numpy, matplotlib, random and scipy libraries. 

This script contains the following functions:
    * pr_theta - converts pseudorapidity values to theta values, where theta
    is the angle to the beam axis.
    * line - a simple function defining the equation of a straight line, this
    is to be used with curve_fit perform straight-line fits.
    * fit_3d - a function which fits straight lines to both the xz and yz
    planes and returns their gradient.

It also contains the following classes:
    * velo - class representing the vertex locator using the left_sensor and
    right_sensor classes, contains methods for calculating hit positions, the
    reconstruction efficiency, the reconstruction efficiency as a function of 
    pseudorapidity and generating particles with uniform randomly sampled
    pseudorapidities and azimuth angles.  '''


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from sensor import left_sensor, right_sensor
from particle import particle
import random
from scipy.optimize import curve_fit



def pr_theta(pseudo_rapidity):
    '''Calculates theta given a pseudo rapidity value
    
    Args:
        pseudo_rapidity - float/ array of psuedo rapidity value(s) to be
        converted.
    
    Returns:
        theta - float/ array of theta values.'''


    theta = 2 * np.arctan(np.exp(-pseudo_rapidity))
    return theta


def line(x, m, c):
    '''Equation of a straight line, to be used with scipy.curve_fit. 
    
    Args:
        x - float/ array, the x values over which the line is defined.
        m - float, the gradient of the line.
        c - float, the y-intercept of the line.
    
    Returns:
        y - float/ array, the y values resulting from given values of m, x and
        c.  '''

    y = m * x + c
    return y
    

def fit_3d(z, x, y, plot = False):
    ''' Fits a straight line to a set of points in 3d space using the curve_fit.

    Args:
        z - list/array, a list or array of z values to be fitted
        x - list/ array, a list or array of x values to be fitted.
        y - list/array, a list or array of y values to be fitted.
        plot - bool, If this is true a plot of the straight line fit will be
        shown.

    Returns: 
        m_x - gradient output by a straight-line fit to data in the xz plane.
        m_y - gradient output by a straight-line fit to data in the yz plane.
        m_x_sigma - error on the value m_x, given by the curve_fit covariance 
        matrix.
        m_y_sigma - error on the value m_y, given by the curve_fit covariance 
        matrix.
    '''

    # Fitting the xz plane
    popt_xz, pcov_xz = curve_fit(line, np.array(z), np.array(x))
    m_x = popt_xz[0]
    c_x = popt_xz[1]
    xz_errors = np.diag(pcov_xz)  # getting errors from covariance matrix
    m_x_sigma = xz_errors[0]
    c_x_sigma = xz_errors[1]
    # theta_x = np.arctan(popt_xz[0])

    # Fitting the xy plane
    popt_yz, pcov_yz = curve_fit(line, z, y)
    m_y = popt_yz[0]
    c_y = popt_yz[1]
    yz_errors = np.diag(pcov_yz)
    m_y_sigma = yz_errors[0]
    c_y_sigma = yz_errors[1]
    # theta_y = np.arctan(popt_yz[0])

    if plot:   # plots each straight line fit
        fig, axs = plt.subplots(2)
        axs[0].plot(z, line(z, *popt_xz), "r--", alpha=0.5, 
            label = f"straight line fit \nm = {popt_xz[0]:.2} $\pm$ {xz_errors[0]:.2}\nc = {popt_xz[1]:.2} $\pm$ {xz_errors[1]:.2}")
        axs[0].scatter(z, x)
        axs[0].set_title("Straight Line Fit")
        axs[0].set_xlabel("z"), axs[0].set_ylabel("x")

        axs[1].plot(z, line(z, *popt_yz), "r--", 
            alpha=0.5, label=f"straight line fit \nm = {popt_yz[0]:.2} $\pm$ {yz_errors[0]:.2} \nc = {popt_yz[1]:.2} $\pm$ {yz_errors[1]:.2}")
        axs[1].scatter(z, y)
        axs[1].set_title("Straight Line Fit")
        axs[1].set_xlabel("z"), axs[1].set_ylabel("y")
        
        axs[0].legend()
        axs[1].legend()

        fig.tight_layout()
    
    
    return m_x, m_y, m_x_sigma, m_y_sigma, c_x, c_y, c_x_sigma, c_y_sigma 
    


class velo:
    ''' A class representing the VELO detector at LHCb.

    Attributes:
        l_sensors - list, a list of left sensor objects, created using the
        left_sensor class
        r_sensors - list, a list of right sensor objects, created using the
        right_sensor class.
        left_sensor_z - list/array, a list or array of z locations of the left
        sensors of the detector.
        right_sensor_z - list/array, a list or array of z locations of the right
        sensors of the detector.
        zs - numpy array, a numpy array containing the z values of both the left 
        and right sensors.  
        particles - array, an array of the generated particle objects.
    '''


    def __init__(self, left_sensor_z, right_sensor_z):
        '''Constructor for the velo class.
        
        Args: 
            left_sensor_z - list/array, a list or array of z locations of the left
            sensors of the detector.
            right_sensor_z - list/array, a list or array of z locations of the right
            sensors of the detector.
        
        Returns:
            None'''


        self.l_sensors = []
        self.r_sensors = []
        self.left_sensor_z = left_sensor_z
        self.right_sensor_z = right_sensor_z
        self.zs = np.concatenate([left_sensor_z, right_sensor_z])

        # creating left sensors at given zs
        for z in self.left_sensor_z:
            l_sensor = left_sensor(z)
            self.l_sensors.append(l_sensor)

        # creating right sensors at given zs
        for z in self.right_sensor_z:
            r_sensor = right_sensor(z)
            self.r_sensors.append(r_sensor)


    def hits(self, particle, hit_efficiency = 98, hit_resolution = 0.012):
        ''' Method for determining the location of sensor hits and the number
        of them. Randomly smears particle position with Gaussian distribution 
        to simulate resolution of detector. 
        
        Args: 
            particle - particle object to be checked for a sensor hit.
            hit_efficiency - float, the percentage chance that the particle is
            detected by the sensor. default: 98.
            hit_resolution - float, the resolution of the sensors in mm.
            default: 0.012 mm.
        
        Returns
            hits - an array of xyz positions of the particle hits.'''


        particle.xy_pos(self.zs)

        hits = []
        

        for i in range(len(self.left_sensor_z)):
            l_x, l_y, l_z = particle.position(self.left_sensor_z[i])
            l_x_smear = np.random.normal(l_x, hit_resolution)  # generating values to smear left x by
            l_y_smear = np.random.normal(l_y, hit_resolution)  # generating values to smear left y by
            l_hits = self.l_sensors[i].in_sensor(l_x, l_y)  # checking if smeared value is in sensor

            # determining if hit
            if (l_hits[0].size > 0 ) & (l_hits[1].size > 0) & (random.choices([True, False], weights = (hit_efficiency, 100 - hit_efficiency))[0]):
                hits.append((l_z, l_x_smear, l_y_smear))
            

        for i in range(len(self.right_sensor_z)):
            r_x, r_y, r_z = particle.position(self.right_sensor_z[i])
            r_x_smear = np.random.normal(r_x, hit_resolution)  # generating values to smear right x by
            r_y_smear = np.random.normal(r_y, hit_resolution)  # generating values to smear right y by 
            r_hits = self.r_sensors[i].in_sensor(r_x, r_y)  # checking if smeared value is in sensor


            if (r_hits[0].size > 0) & (random.choices([True, False], weights = (hit_efficiency, 100 - hit_efficiency))[0]):
                hits.append((r_z, r_x_smear, r_y_smear))


        return hits
    

    def uniform_particle_generate(self, phi_range: tuple, eta_range: tuple, number: int, start_pos_range=(0,0,0)):
        '''Generates a given number of particle objects with pseudo rapidity in the
        given phi and eta ranges.

        Args:
            phi_range - tuple of floats, the range of phis that particles can be
            generated with.
            eta_range - tuple of floats, the range of etas that particles can be
            generated with.
            number - int, the number of particle objects to be created.
            pseudo_rapidity - int/ float, the pseudo rapidity of the desired particles.
            start_pos_range - tuple of lists, the x, y and z ranges in which
            random uniform values for the start position are generated.
            default: (0, 0, 0).

        Returns:
            self.particles - list, A list of particle objects, created using the 
            particle class.
            phi_values - numpy array, An array of the phi values of the generated
            particles.
            eta_values - numpy array, An array of the eta values of the generated 
            particles.
        '''


        phi_values = np.random.uniform(*phi_range, size=number)
        eta_values = np.random.uniform(*eta_range, size=number)
        start_pos_vals = np.repeat([[0, 0, 0]], repeats=number, axis=0)  # creating a 2d array of start pos

        if start_pos_range != (0,0,0):  # generates random start pos values
            x_start = np.random.uniform(start_pos_range[0][0], start_pos_range[0][1], size=number)
            y_start = np.random.uniform(start_pos_range[1][0], start_pos_range[1][1], size=number)
            z_start = np.random.uniform(start_pos_range[2][0], start_pos_range[2][1], size=number)
            start_pos_vals = np.array(np.transpose([x_start, y_start, z_start]))  # transposing our vector for ease
            
        self.particles = []

        for count, eta in enumerate(eta_values):
            theta = pr_theta(eta)  # converts pseudorapidity to angle to beam
            start = start_pos_vals[count]
            ys = self.zs * np.tan(theta) 
            xs = self.zs * ((np.tan(theta))/(np.tan(phi_values[count])))
            part_obj = particle((xs[0] ,ys[0], self.zs[0]), start)  # create a particle with given values
            part_obj.xy_pos(self.zs)
            self.particles.append(part_obj)

        return self.particles, phi_values, eta_values


    def recon_eff(self):
        ''' Calculates the reconstruction efficiency for a given set of particles.

        Args:
            None.

        Returns:
            efficiency - float, the reconstruction efficency of the given set 
            of particles.
            reconstructed_hits - array, an array containing arrays of the hit
            positions for each reconstructed particle.
        '''


        reconned = 0
        reconstructed_hits = []
        for par in self.particles:
            hts = self.hits(par)

            # checking if a particle is reconstructed
            if len(hts) >= 3:
                reconned += 1
                reconstructed_hits.append(hts)

        efficiency = reconned / len(self.particles)
        return efficiency, reconstructed_hits


    def recon_eta(self, pr_range, start_pos):
        ''' Calculates the reconstruction efficiency of a given velo object.

        Args:
            pr_range - tuple of floats, the range which psuedorapidity values
            can lie within,.
            start_pos - tuple of floats, the position at which the the particle
            is generated. 
        
        Returns:
            etas - array of floats, the pseudorapidity values at which the
            particles are generated.
            effs - array of floats, the reconstruction efficiencies that each
            set of particles has.
        '''


        etas = np.linspace(*pr_range, 100)  # generating uniformly spaced etas
        effs = []
        for eta in etas:
            particles = self.uniform_particle_generate((0, 2 * np.pi), (eta, eta), 100, start_pos)[0]  # generating particles
            eff, recon_hits = self.recon_eff()
            effs.append(eff)
            
        return etas, effs


if __name__ == "__main__":
    pass
