# Year 3 Computing Labs Project - Modelling VELO at the LHCb
## Solomon Sanderson

This repository contains the code from my Year 3 Particle Physics Computing laboratory. 

In this laboratory I aimed to model the vertex locator (VELO) at LHCb, CERN. This can be seen on the diagram below:

<img src="media\cern detector.png">

This simulation accounted for numerous effects:
* The resolution of the detector, simulated using random gaussian numbers. 
* Hit efficiency - simulated using random numbers.
* Detector Geometry. 

For more detail please see my report in this [PDF](https://drive.google.com/file/d/1sTPed9F4mVLHar0WHeR5y1KZ0L9g6Hjo/view?usp=sharing).

The files are as follows:
* `main.py` - This is where the code should be run from, this imports all of the other modules.

* `particle.py` - Defines the particle class, this holds the particles position and momentum values, it also holds methods for setting the momentum of the particle. It also holds methods which let us get the transverse momenta, impact parameter and their respective resolutions.

* `cube_plot.py` - Plots cubes for making the 3-dimensional plots of the vertex locator detectors sensors, seen in red and blue on the plots. 

* `sensor.py` - Contains classes for the left and right sensors of the VELO detector. Contains methods which return whether or not a particle is in a sensor. Also contains a function for creating the 3d plots, some examples are below: 

<img src="media\tracks_varied_eta.png">
<img src="media\tracks_varied_phi.png">



* `velo.py` - Contains methods for generating particles with given angles and momenta, calculating the reconstruction efficiency and fitting lines in 3d to determine particle paths.