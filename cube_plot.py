''' cube_plot.py

This script creates matplotlib cube objects which are used in the plotting of
my 3d sensor model.
https://stackoverflow.com/questions/49277753/python-matplotlib-plotting-cuboids

This code contains 2 functions: 
    * cuboid_data2 - returns a 3 dimensional array which defines each of the
      vertices of the cuboids

'''


from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import matplotlib.pyplot as plt


def cuboid_data2(o, size=(1,1,1)):
    '''
    Creates a 3 dimensional array which defines the vertices of the cuboids. 

    Args: 
        o - tuple/list, The positions of the cuboids, passed in the form 
        (x, y, z).
        size - tuple/list, The size of the cuboids to be generated defaults to
        (1,1,1).
    
    Returns:
        X - 2d array, contains points defining the 4 corners of each face of
        cuboid.

    '''
    # 3d array/ ,atric representing the corners of the cube
    X = [[[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0]],
         [[0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 0, 0]],
         [[1, 0, 1], [1, 0, 0], [1, 1, 0], [1, 1, 1]],
         [[0, 0, 1], [0, 0, 0], [0, 1, 0], [0, 1, 1]],
         [[0, 1, 0], [0, 1, 1], [1, 1, 1], [1, 1, 0]],
         [[0, 1, 1], [0, 0, 1], [1, 0, 1], [1, 1, 1]]]
    X = np.array(X).astype(float)  # changing vakues to floats
    for i in range(3):
        X[:,:,i] *= size[i]
    X += np.array(o)  # Adding position to X
    return X


def plotCubeAt2(positions,sizes=None,colors=None, **kwargs):
    '''
    Creates a matplotlib cuboid object of the given size, position and color
    (if given). 

    Args:
        positions - 3d array, the position at which to create the cuboid.
        sizes - 3d array, the sizes of the cuboids that will be generated.
        colors - array of strings, the colours of each of the cuboid objects
    
    Returns:
        Poly3DCollection - Collection of matplotlib 3D objects.
    '''
    if not isinstance(colors,(list,np.ndarray)): colors=["C0"]*len(positions)  # checks colors is a tuple of a list and an array.
    if not isinstance(sizes,(list,np.ndarray)): sizes=[(1,1,1)]*len(positions)  # as above for sizes.
    g = []
    for p,s,c in zip(positions,sizes,colors):  # loops over positions, sizes and colours
        g.append( cuboid_data2(p, size=s)) 
    return Poly3DCollection(np.concatenate(g),  
                            facecolors=np.repeat(colors,6), **kwargs)

                            
if __name__ == "__main__":   
    pass