''' particle.py

Allows the create particle objects to be used in the modelling of the sensor.

This script requires the numpy library. 

It contains the following classes:
    * particle - creates a particle object which has a given momentum and
    start position. The momentum can be set to the desired magnitude. Has 
    methods to determine the transverse momentum, its error and the impact
    parameter.
'''

import numpy as np


class particle:
    '''A class for a particle to be detected by the sensors
    
    Attributes:
        p_arr - array, an array of the particles initial momentum values.
        p_x - float, the initial momentum in the x direction.
        p_y - float, the initial momentum in the y direction.
        p_z - float, the initial momentum in the z direction.
        x_0 - float, the initial position in the x direction.
        y_0 - float, the initial position in the y direction.
        z_0 - float, the initial position in the z direction.
        '''


    def __init__(self, momentum, start):
        '''Constructor for particle class
        
        Args:
            momentum: tuple/list - a momentum vector of the form (p_x, p_y, p_z)
            start: tuple/list - the starting position of the particle that we measure (x, y, z))
        
        Returns:
            None
        '''

        # (start)
        self.p_arr = np.array(momentum)
        self.p_x, self.p_y, self.p_z = momentum
        self.x_0, self.y_0, self.z_0 = start


    def xy_pos(self, z_arr):
        '''Takes a z_arr to calculate the nearest z point to the starting 
        position. Uses this to calculate the trajectory of the particle and the
        x, y, z values of its path.
        
        Args:
            z_arr - numpy array, array of z positions of sensors.
            
        Returns:
            None.'''

        # print(self.x_0)
        
        self.x_mult = self.p_x/self.p_z  # gradient 
        self.y_mult = self.p_y/self.p_z  # gradient
        self.z_arr = z_arr 
        

        # self.x_arr = (z_arr * self.x_mult) + self.x_0
        # self.y_arr = (z_arr * self.y_mult) + self.y_0
        # shifted_z = z_arr + self.z_0

        # self.delta_x_z = (self.x_mult * self.z_0) + self.x_0
        # self.delta_y_z = (self.y_mult * self.z_0) + self.y_0
        
        # self.x_arr = (z_arr * self.x_mult)
        # self.y_arr = (z_arr * self.y_mult)
        self.x_arr = ((z_arr - self.z_0) * self.x_mult) + self.x_0 #+ self.delta_x_z
        self.y_arr = ((z_arr - self.z_0) * self.y_mult) + self.y_0 #+ self.delta_y_z
        # print(self.x_arr[0], self.y_arr[0], z_arr[0])
        # self.z_0_x = -(self.x_0 + self.delta_x_z)/self.x_mult
        # self.z_0_y = -(self.y_0 + self.delta_y_z)/self.y_mult
        # # print(f"z_0 {self.z_0_x, self.z_0_y}")
        # self.x_0 = self.x_0 + delta_x_z
        # self.y_0 = self.y_0 + delta_y_z
        


    def position(self, z_value):
        ''' Returns the x and y values of a particles position for a given 
        z-value.
        
        Args:
            z_value - float/ int, the z value at which the position of the
            particle should be calculated.
            
        Returns:
            (x, y, z) - tuple, a tuple of position values at the specified
            z position.'''


        index = np.where(z_value == self.z_arr)
        return self.x_arr[index], self.y_arr[index], self.z_arr[index] 

        
    def set_pmag(self, momentum_magnitude):
        '''
        Sets the magnitude of the particle to the given value.

        Args:
            momentum_magnitude - float, the desired magnitude of the particle
            momentum.
        Returns:
            None.
        '''


        current_pmag = np.linalg.norm(self.p_arr)
        unit_vector = (1 / current_pmag) * self.p_arr
        new_p = unit_vector * momentum_magnitude
        self.p_arr = np.array(new_p)  # consider removing the individual pxpypz part 
        self.xy_pos(self.z_arr)
        # self.p_x, self.p_x, self.p_z = self.p_arr


    def transverse_momentum(self, m_x, m_y, p_tot):
        '''
        Calculates the transverse momentum of the particle.

        Args:
            m_x - float, the gradient of the line in the xz plane fitted by
            scipy.curvefit.
            m_y - float, the gradient of the line in the yz plane 
            fitted by scipy.curvefit.
            p_tot - float, the total momentum of the particle. 
        Returns
            p_t - float, the transverse momentum of the particle.The magnitude
            of the momentum in the xy plane.
        '''


        m_squared = m_x ** 2 + m_y ** 2
        p_t = p_tot * np.sqrt(m_squared/(m_squared + 1))
        return p_t


    def p_resolution(self, m_x, m_y, m_x_error, m_y_error):
        '''
        Calculates the resolution of the transverse momentum values.

        Args:
            m_x - float, gradient of our straight-line fit in the xz plane.
            m_y - float, gradient of our straight line fit in the yx plane.
            m_x_error - float, the error on m_x given by the covariance matrix
            from scipy.curve_fit.
            m_y_error - float, the error on m_y given by the covariance matrix
            from scipy.curve_fit.

        Returns:
            sigma_p_t - float, the error on the transverse momentum of the
            particle.
        '''

        # print(m_x, m_y, m_x_eror, m_y_error)
        # sigma_p_t = (2 * np.sqrt((m_x_error / m_x) ** 2 + (m_y_error / m_y) ** 2))/(1 + 2 * np.sqrt((m_x_error / m_x) ** 2 + (m_y_error / m_y) ** 2))
        p_x = m_x * self.p_arr[2]
        p_y = m_y * self.p_arr[2]
        # print(self.p_arr)
        p_tot = (p_x ** 2 + p_y ** 2 + self.p_arr[2] ** 2) ** 0.5


        # print(self.p_arr[0])
        # print(f"p_x{p_x}")
        # print(f"p_y{p_y}, p_tot{p_tot}")
        sigma_p_x =  m_y_error * np.abs(self.p_arr[2])
        sigma_p_y =  m_x_error * np.abs(self.p_arr[2])
        # print(f"self.p_arr[2] {self.p_arr[2]}")

        sigma_p_t = 0.5 * ( 2 *((sigma_p_x / p_x) + (sigma_p_y / p_y)))/((p_x ** 2 + p_y ** 2) ** 0.5)
        # print(sigma_p_t)
        return sigma_p_t


    def impact_parameter(self, m_x, m_y, c_x, c_y):
        '''
        Calculates the impact parameter (closest distance of approach to the 
        true origin) of our particle.
        
        Args:
            None.
        
        Returns:
            dist - float, the shortest distance between the particles true origin and
            the particle track.
            '''
        
        self.x_fitted = m_x * self.z_arr + c_x
        self.y_fitted = m_y * self.z_arr + c_y

        x_1, x_2 = self.x_fitted[0], self.x_fitted[1]
        y_1, y_2 = self.y_fitted[0], self.y_fitted[1]
        z_1, z_2 = self.z_arr[0], self.z_arr[1]

        t = -((x_2 - self.x_0) * (x_1 - x_2) + (y_2 - self.y_0) * (y_1 - y_2) + (z_2 - self.z_0) * (z_1 - z_2))/((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2)
        b_x = t * (x_1 - x_2) + (x_2 - self.x_0)
        b_y = t * (y_1 - y_2) + (y_2 - self.y_0)
        b_z = t * (z_1 - z_2) + (z_2 - self.z_0)
        dist = np.linalg.norm([b_x, b_y, b_z])
        return dist

    
    def ip_resolution(self, m_x, m_y, m_x_sigma, m_y_sigma, c_x, c_y, c_x_sigma, c_y_sigma ):
        '''
        Calculates the resolution of the impact parameter.

        Args:
            m_x - float, the gradient returned by the straight line fit in the 
            xz plane.
            m_y - float, the gradient returned by the straight line fit in the 
            yz plane.
            m_x_sigma - float, the error on the gradient from the covariance
            matrix returned by the straight line fit.
            m_y_sigma - float, the error on the gradient from the covariance
            matrix
        '''
        sigma_fit_x = self.z_arr * (m_x_sigma / m_x) + c_x_sigma
        sigma_fit_y = self.z_arr * (m_y_sigma / m_y) + c_y_sigma

        x_1, x_2 = self.x_fitted[0], self.x_fitted[1]
        y_1, y_2 = self.y_fitted[0], self.y_fitted[1]
        z_1, z_2 = self.z_arr[0], self.z_arr[1]
        t = -((x_2 - self.x_0) * (x_1 - x_2) + (y_2 - self.y_0) * (y_1 - y_2) + (z_2 - self.z_0) * (z_1 - z_2))/((x_1 - x_2) ** 2 + (y_1 - y_2) ** 2 + (z_1 - z_2) ** 2)
        b_x = t * (x_1 - x_2) + (x_2 - self.x_0)
        b_y = t * (y_1 - y_2) + (y_2 - self.y_0)
        b_z = t * (z_1 - z_2) + (z_2 - self.z_0)

        sigma_t = (sigma_fit_x ** 2 + sigma_fit_y ** 2)/( 2 * ((2 * sigma_fit_x)/(x_1 - x_2) + (2 * sigma_fit_y)/(y_1 - y_2)))
        sigma_b_x = (sigma_t)/(t) + (2 * sigma_fit_x)/(x_1 - x_2) + sigma_fit_x
        sigma_b_y = (sigma_t)/(t) + (2 * sigma_fit_y)/(y_1 - y_2) + sigma_fit_y
        sigma_xyz = (2 * sigma_b_x)/(b_x) + (2 * sigma_b_y)/(b_y)
        sigma_dist = 0.5  * ((sigma_xyz)/(b_x ** 2 + b_y ** 2 + b_z ** 2))
        return sigma_dist   
    

if __name__ == "__main__":
    pass