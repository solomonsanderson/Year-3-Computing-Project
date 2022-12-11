'''Creating a particle class'''


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


        self.p_arr = np.array(momentum)
        self.p_x, self.p_y, self.p_z = momentum
        self.x_0, self.y_0, self.z_0 = start


    def calibrate(self, z_arr):
        '''Takes a z_arr to calculate the nearest z point to the starting 
        position. Uses this to calculate the trajectory of the particle and the
        x, y, z values of its path.
        
        Args:
            z_arr - numpy array, array of z positions of sensors.
            
        Returns:
            None.'''


        nearest_z = z_arr[(np.abs(z_arr - self.z_0)).argmin()]
        z_diff = nearest_z - self.z_0  # distance between starting point and first sensor
        
        self.x_mult = self.p_x/self.p_z
        self.y_mult = self.p_y/self.p_z

        self.z_arr = z_arr
        self.x_arr = z_arr * self.x_mult
        self.y_arr = z_arr * self.y_mult 
        

    def position(self, z_value):
        '''Calculates the position of the particle in the x and y plane and
        returns the coordinates in the format (x, y, z).
        
        Args:
            z_value - float/ int, the z value at which the position of the
            particle should be calculated.
            
        Returns:
            (x, y, z) - tuple, a tuple of position values at the specified
            z position.'''


        # print(z_value)
        index = np.where(z_value == self.z_arr)
        # print(self.x_arr[index], self.y_arr[index], self.z_arr[index])
        return self.x_arr[index], self.y_arr[index], self.z_arr[index] 

        
    def set_pmag(self, momentum_magnitude):
        '''
        
        '''


        current_pmag = np.linalg.norm(self.p_arr)
        unit_vector = (1 / current_pmag) * self.p_arr
        new_p = unit_vector * momentum_magnitude
        # print(new_p)
        # print(np.linalg.norm(new_p))
        self.p_arr = np.array(new_p)  # consider removing the individual pxpypz part 
    

    def transverse_momentum(self, m_x, m_y, p_tot):
        '''
        
        '''


        # px_over_py = np.tan(theta_x)
        # py_over_pz = np.tan(theta_y)
        m_squared = m_x ** 2 + m_y ** 2
        p_t = p_tot * np.sqrt(m_squared/(m_squared + 1))
        return p_t



    # def p_resolution(self, m_x, sigma_mx, m_y, sigma_my):
    #     theta_x, theta_y = np.arctan(m_x), np.arctan(m_y)
    #     p_x = p_xz * np.sin(theta_x)
    #     p_y = p_yz * np.sin(theta_y)
    #     p_t = np.linalg(np.array([p_x, p_y]))
    #     sigma_pt = (p_x *  sigma_mx *(np.cos(theta_x)) ** 3) + (p_y * sigma_my * (np.cos(theta_y)) ** 3)

    def p_resolution(self, m_x, m_y, m_x_error, m_y_error):
        '''
        '''

        # print(m_x, m_y, m_x_error, m_y_error)
        sigma_p_t = (2 * np.sqrt((m_x_error / m_x) ** 2 + (m_y_error / m_y) ** 2))/(1 + 2 * np.sqrt((m_x_error / m_x) ** 2 + (m_y_error / m_y) ** 2))
        # print(sigma_p_t)
        return sigma_p_t



if __name__ == "__main__":
    pass