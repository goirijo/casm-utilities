import casmutils as cu
import numpy as np


def _conventional_rocksalt_transformation_matrix():
    """Return
    0 1 1
    1 0 1
    1 1 0
    Returns
    -------
    np.array

    """
    m = np.ones([3, 3])
    for i in range(3):
        m[i, i] = 0

    return m

def _final_lattice(block_dims, shift_vals):
    """Creates the final lattice of the Wadsley-Roth
    structure. 

    Parameters
    ----------
    block_dims : pair of int, dimensions of block unit in terms of octahedra
    shift_vals : pair of int, define the lateral and vertical shear

    Returns
    -------
    np.array

    """
    bx,by=block_dims
    sx,sy=shift_vals
    S = np.array([[2 * bx - 1, 2 * sy, 0], [2 * sx, 2 * by - 1, 0],
                  [1, 1, 2]])
    return S

def single_block_wadsley_roth(block_dims, shift_vals, species, nearest_neighbor_distance):
    """Generates a Wadsley-Roth structure where there
    is only a single tiling unit (e.g. 3x3 octahedra)

    Parameters
    ----------
    block_dims : pair of int, dimensions of block unit in terms of octahedra
    shift_vals : pair of int, define the lateral and vertical shear
    species : pair of str, octahedral center goes first
    nearest_neighbor_distance : float

    Returns
    -------
    xtal.Structure

    """
    m=_conventional_rocksalt_transformation_matrix()
    S=_final_lattice(block_dims,shift_vals)
    trans = np.dot(np.linalg.inv(m), S)

    nn = nearest_neighbor_distance
    dnn = 2 * nn
    rs = cu.xtal.RockSaltToggler.relative_to_primitive(trans, species, 0, nn)

    for j in range(block_dims[0]):
        for k in range(block_dims[1]):
            coord = cu.xtal.Coordinate([dnn * j, dnn * k, 0])
            rs.activate(coord)

    return rs.structure()

def main():
    block_dims=[6,4]
    shift_vals=[2,0]

    tilel,tilew=block_dims
    mshift,nshift=shift_vals

    rs=single_block_wadsley_roth(block_dims,shift_vals,["Nb","O"],2.0)
    rs.to_poscar("shear" + str(tilel) + "x" + str(tilew) + "_" + str(nshift) + "_" +
                 str(mshift) + ".vasp")
    

if __name__=="__main__":
    main()
