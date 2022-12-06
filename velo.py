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



def pr_theta(pseudo_rapidity):
    '''Calculates theta given a pseudo rapidity value
    
    Args:
        pseudo_rapidity - float/ array of psuedo rapidity value(s) to be
        converted.
    
    Returns:
        theta - float/ array of theta values.'''


    theta = 2 * np.arctan(np.exp(-pseudo_rapidity))
    return theta


def particle_generate(number, pseudo_rapidity, z):
    '''Generates a given number of particle objects with pseudo rapidity at 
    given z values.
    
    Args:
        number - int, the number of particle objects to be created.
        pseudo_rapidity - int/ float, the pseudo rapidity of the desired particles.
        z - array/ int, the z positions of the sensors.
    
    Returns:
        particles - list, A list of particle objects, created using the 
        particle class.
    '''


    theta = pr_theta(pseudo_rapidity)
    ys = z * np.tan(theta)

    particles = []

    
    while len(particles) < number:
        phi = np.random.uniform(0, 2 * np.pi)
        xs = z * ((np.tan(theta))/(np.tan(phi)))
        part_obj = particle((xs[0] ,ys[0], z[0]), (0,0,0))
        part_obj.calibrate(z)
        particles.append(part_obj)
    
    return particles


def uniform_particle_generate(phi_range, eta_range, number, z):
    '''
    
    '''


    phi_values = np.linspace(phi_range[0], phi_range[1], number)
    eta_values = np.linspace(eta_range[0], eta_range[1], number)
    particles = []


    for phi, eta in list(zip(phi_values, eta_values)):
        theta = pr_theta(eta)
        ys = z * np.tan(theta)
        xs = z * ((np.tan(theta))/(np.tan(phi)))
        part_obj = particle((xs[0] ,ys[0], z[0]), (0,0,0))
        part_obj.calibrate(z)
        particles.append(part_obj)

    
    # print(particles)
    return particles


def recon_eff(particles, velo_obj):
    ''' Calculates the reconstruction efficiency for a given set of particles.
    
    Args:
        particles - list, A list of particle objects.

    Returns:
        efficiency - float, the reconstruction efficency of the given set 
        of particles.
    '''


    reconned = 0
    reconstructed_hits = []
    for par in particles:
        hits, number = velo_obj.hits(par)
        if number >= 3:
            reconned += 1
            reconstructed_hits.append(hits)

        # ax3d.plot(par.z_arr, par.x_arr,  par.y_arr,  marker=None)

    efficiency = reconned / len(particles)
    return efficiency, reconstructed_hits

def line(x, m, c):
    y = m * x + c
    return y
    
def fit_3d(z, x, y):
    ''' Fits a straight line to a set of points in 3d space using the First 
    principle component.

    Args:
        x - list/ array, a list or array of x values to be fitted.
        y - list/array, a list or array of y values to be fitted.
        z - list/array, a list or array of z values to be fitted
    
    returns: 
        linepts - numpy array, an array of points for a straight line
        which has been fit to the data.
    '''

    print(z,x)
    fig, axs = plt.subplots(2)
    popt_xz, pcov_xz = curve_fit(line, np.array(z), np.array(x))
    xz_errors = np.diag(pcov_xz)
    axs[0].plot(z, line(z, *popt_xz), "r--", alpha=0.5, label = f"straight line fit \nm = {popt_xz[0]:.2} $\pm$ {xz_errors[0]:.2}\nc = {popt_xz[1]:.2} $\pm$ {xz_errors[1]:.2}")
    axs[0].scatter(z, x)
    axs[0].set_title("Straight Line Fit")
    axs[0].set_xlabel("z"), axs[0].set_ylabel("x")
    zx_angle = np.arctan(popt_xz[1])

    popt_yz, pcov_yz = curve_fit(line, z, y)
    yz_errors = np.diag(pcov_yz)
    axs[1].plot(z, line(z, *popt_yz), "r--", alpha=0.5, label=f"straight line fit \nm = {popt_yz[0]:.2} $\pm$ {yz_errors[0]:.2} \nc = {popt_yz[1]:.2} $\pm$ {yz_errors[1]:.2}")
    axs[1].scatter(z, y)
    axs[1].set_title("Straight Line Fit")
    axs[1].set_xlabel("z"), axs[1].set_ylabel("y")
    zy_angle = np.arctan(popt_yz[1])
    print(np.degrees(zx_angle), np.degrees(zy_angle))
    axs[0].legend()
    axs[1].legend()
    fig.tight_layout()
    


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


        for z in self.left_sensor_z:
            l_sensor = left_sensor(z)
            self.l_sensors.append(l_sensor)

        for z in self.right_sensor_z:
            r_sensor = right_sensor(z)
            self.r_sensors.append(r_sensor)


    def hits(self, particle, hit_efficiency = 98, hit_resolution = 0.012):
        ''' Method for determining the location of sensor hits and the number
        of them. Randomly smears particle position with Gaussian distribution 
        to simulate resolution of detector. 
        
        Args: particle particle object to be checked for a sensor hit.
            hit_efficiency - float, the percentage chance that the particle is
            detected by the sensor. default: 98.
            hit_resolution - float, the resolution of the sensors. '''


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
    pass
