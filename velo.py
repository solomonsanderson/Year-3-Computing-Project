'''creating a class for the velo sensor'''


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
        c - float, the y-intercept of the line '''

    y = m * x + c
    return y
    

def fit_3d(z, x, y, plot = False):
    ''' Fits a straight line to a set of points in 3d space using the First 
    principle component.

    Args:
        z - list/array, a list or array of z values to be fitted
        x - list/ array, a list or array of x values to be fitted.
        y - list/array, a list or array of y values to be fitted.
    
    returns: 
        linepts - numpy array, an array of points for a straight line
        which has been fit to the data.
    '''


    popt_xz, pcov_xz = curve_fit(line, np.array(z), np.array(x))
    m_x = popt_xz[0]
    xz_errors = np.diag(pcov_xz)
    m_x_sigma = xz_errors[0]
    # theta_x = np.arctan(popt_xz[0])

    popt_yz, pcov_yz = curve_fit(line, z, y)
    m_y = popt_yz[0]
    # theta_y = np.arctan(popt_yz[0])
    yz_errors = np.diag(pcov_yz)
    m_y_sigma = yz_errors[0]

    if plot:
        fig, axs = plt.subplots(2)
        axs[0].plot(z, line(z, *popt_xz), "r--", alpha=0.5, label = f"straight line fit \nm = {popt_xz[0]:.2} $\pm$ {xz_errors[0]:.2}\nc = {popt_xz[1]:.2} $\pm$ {xz_errors[1]:.2}")
        axs[0].scatter(z, x)
        axs[0].set_title("Straight Line Fit")
        axs[0].set_xlabel("z"), axs[0].set_ylabel("x")

        axs[1].plot(z, line(z, *popt_yz), "r--", alpha=0.5, label=f"straight line fit \nm = {popt_yz[0]:.2} $\pm$ {yz_errors[0]:.2} \nc = {popt_yz[1]:.2} $\pm$ {yz_errors[1]:.2}")
        axs[1].scatter(z, y)
        axs[1].set_title("Straight Line Fit")
        axs[1].set_xlabel("z"), axs[1].set_ylabel("y")
        
        axs[0].legend()
        axs[1].legend()

        fig.tight_layout()
    
    # print(np.degrees(zx_angle), np.degrees(zy_angle))
    
    return m_x, m_y, m_x_sigma, m_y_sigma
    


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

        particle.xy_pos(self.zs)
        
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
            z - array/ int, the z positions of the sensors.

        Returns:
            particles - list, A list of particle objects, created using the 
            particle class.
            phi_values - numpy array, An array of the phi values of the generated
            particles.
            eta_values - numpy array, An array of the eta values of the generated 
            particles
        '''


        phi_values = np.random.uniform(*phi_range, size=number)
        eta_values = np.random.uniform(*eta_range, size=number)
        # print(eta_values)

        start_pos_vals = np.repeat([[0, 0, 0]], repeats=number, axis=0)


        if start_pos_range != (0,0,0):
            x_start = np.random.uniform(start_pos_range[0][0], start_pos_range[0][1], size=number)
            y_start = np.random.uniform(start_pos_range[1][0], start_pos_range[1][1], size=number)
            z_start = np.random.uniform(start_pos_range[2][0], start_pos_range[2][1], size=number)
            start_pos_vals = np.array(np.transpose([x_start, y_start, z_start]))
            # print(f"starpv{start_pos_vals}")
        # print(start_pos_vals)
        
        # print(f"phis = {phi_values}")
        # print(f"etas = {eta_values}")
        self.particles = []
        # print(start_pos_vals)

        for count, eta in enumerate(eta_values):
            # print(start_pos_vals[count])
            # try enumerate here with x,y,z
            # need to only give one value for x y z in each iteration
            # print(f"phi{phi}, eta{eta}")
            theta = pr_theta(eta_values[count])
            # print(f"start{start}")
            start = start_pos_vals[count]
            ys = self.zs * np.tan(theta)
            xs = self.zs * ((np.tan(theta))/(np.tan(phi_values[count])))
            part_obj = particle((xs[0] ,ys[0], self.zs[0]), start)
            part_obj.xy_pos(self.zs)
            self.particles.append(part_obj)


        return self.particles, phi_values, eta_values


    def recon_eff(self):
        ''' Calculates the reconstruction efficiency for a given set of particles.

        Args:
            particles - list, A list of particle objects.

        Returns:
            efficiency - float, the reconstruction efficency of the given set 
            of particles.
        '''


        reconned = 0
        reconstructed_hits = []
        for par in self.particles:
            hts = self.hits(par)
            # if hts != []:
            #     print(hts)

            if len(hts) >= 3:
                reconned += 1
                reconstructed_hits.append(hts)

            # ax3d.plot(par.z_arr, par.x_arr,  par.y_arr,  marker=None)
        efficiency = reconned / len(self.particles)
        return efficiency, reconstructed_hits


    def recon_eta(self, pr_range, start_pos):
        ''' Calculates the reconstruction efficiency of a given velo object 
        '''
        etas = np.linspace(*pr_range, 100)
        effs = []
        for eta in etas:
            # particles = particle_generate(10, eta, z)
            # print(eta)
            # print(start_pos)
            particles = self.uniform_particle_generate((0, 2 * np.pi), (eta, eta), 10, start_pos)[0]
            eff, recon_hits = self.recon_eff()
            effs.append(eff)
            
        # print(etas, effs)
        return etas, effs

if __name__ == "__main__":
    pass
