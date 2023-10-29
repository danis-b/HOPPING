import numpy as np
import json
from io import StringIO
from typing import OrderedDict
from datetime import datetime
import argparse

#https://triqs.github.io/tprf/latest/reference/python_reference.html#wannier90-tight-binding-parsers
def parse_hopping_from_wannier90_hr_dat(filename):

    with open(filename, 'r') as fd:
        lines = fd.readlines()

    lines.pop(0) # pop time header

    num_wann = int(lines.pop(0))
    nrpts = int(lines.pop(0))

    nlines = int(np.ceil(float(nrpts / 15.)))

    deg = np.array([])
    for line in lines[:nlines]:
        deg = np.concatenate((deg, np.loadtxt(StringIO(line), dtype=int, ndmin=1)))
    
    assert( deg.shape == (nrpts,) )

    hopp = "".join(lines[nlines:])
    hopp = np.loadtxt(StringIO(hopp))

    assert( hopp.shape == (num_wann**2 * nrpts, 7) )
    
    R = np.array(hopp[:, :3], dtype=int) # Lattice coordinates in multiples of lattice vectors
    nm = np.array(hopp[:, 3:5], dtype=int) - 1 # orbital index pairs, wannier90 counts from 1, fix by remove 1
    
    n_min, n_max = R.min(0), R.max(0)

    t_re = hopp[:, 5]
    t_im = hopp[:, 6]
    t = t_re + 1.j * t_im # complex hopping amplitudes for each R, mn (H(R)_{mn})

    # -- Dict with hopping matrices

    r_dict = OrderedDict()
    hopp_dict = {}
    for idx in range(R.shape[0]):
        r = tuple(R[idx])

        if r not in r_dict:
            r_dict[r] = 1
        else:
            r_dict[r] += 1

        if r not in hopp_dict:
            hopp_dict[r] = np.zeros((num_wann, num_wann), dtype=complex)

        n, m = nm[idx]
        hopp_dict[r][n, m] = t[idx]

    # -- Account for degeneracy of the Wigner-Seitz points
    
    for r, weight in zip(list(r_dict.keys()), deg):
        hopp_dict[r] /= weight

    return hopp_dict, num_wann, n_min, n_max



def coordination_sort(atom, num_atoms, n_min, n_max, cell_vectors, positions):

    n_size = n_max - n_min + 1  # Plus 1 for 0th
    num_points = num_atoms * n_size.prod()

    # Index filling
    index_ = np.stack(np.unravel_index(np.arange(num_points),
                                       (*n_size, num_atoms)),
                      axis=-1)
    index_ += [*n_min, 0]

    r = index_[:, :3] @ cell_vectors + positions[index_[:, 3]] - positions[atom]
    radius_ = (r * r).sum(-1)**0.5

    # Sort radius descending
    idx = radius_.argsort()
    return radius_[idx], index_[idx]


def main():
    """
    This script converts input_hr.dat hamiltonian from the wannier90 package to block-divided 
    form between orbitals and prints these blocks depending on the geometry given by input.json file.

    Usage: python HOPPING.py input out.dat

    input.json contains the following information:
    "cell_vectors": (3x3)(dfloat) array of unit cell vectors (in Ang);
    "number_of_atoms": (int) total number of atoms, parametrized by Wannier functions;
    "max_sphere_num": (int) maximum number of coordination sphere to print the hamiltonian blocks between pair of atoms;
    "print_complex": (0 or 1) print the imaginary part of hamiltonian blocks;
    "positions": (3 x number_of_atoms)(dfloat) array of cartesian coordinates of these atoms (in Ang);
    "wanniers":  (number_of_magnetic_atoms) (int) array of Wannier functions numbers of these atoms;
    """
    
    parser = argparse.ArgumentParser(prog='HOPPING.py', usage='%(prog)s input (name of input_hr.dat and input.json) output')
    parser.add_argument("input_file_name")
    parser.add_argument("output_file_name")
    args = parser.parse_args()


    print("Program HOPPING.x v.3.0 (python) starts on ", datetime.now())
    print('=' * 69)

    hops, num_wannier_funcs, n_min, n_max = parse_hopping_from_wannier90_hr_dat(f'{args.input_file_name}_hr.dat') 
    n_size = n_max - n_min + 1  # Plus 1 for 0th

    Ham_R = np.zeros((n_size[0], n_size[1], n_size[2], num_wannier_funcs, num_wannier_funcs), dtype='c16')
    
    for r in hops.keys():
        r_idx = np.array(r)
        Ham_idx = hops.get(r)
        
        for m in range(num_wannier_funcs):
            for n in range(num_wannier_funcs):
                Ham_R[r_idx[0] + n_max[0], r_idx[1] + n_max[1], r_idx[2] + n_max[2], m, n] = Ham_idx[m,n]

    

    # =======================================================================
    # Read information from input file in json format
    with open(f'{args.input_file_name}.json') as fp:
        data = json.load(fp)

    num_atoms = data['number_of_atoms']  # M
    positions = np.array(data['positions'])  # array[M, 3]
    wanniers = np.array(data['wanniers'])  # array[M]
    cell_vectors = np.array(data['cell_vectors'])  # array[3, 3]
    max_sphere_num = data['max_sphere_num']
    print_complex = data['print_complex']  # Pretty repr

    # Some checks
    assert num_atoms == len(wanniers) == len(positions)
    assert num_wannier_funcs == wanniers.sum()

    # ===============================================================================
    # Print data to stdout and out-file

    print('Crystal structure of system:')
    print('crystal axes: (cart. coord. )')

    for cell in cell_vectors:
        print('{:.6f} {:.6f} {:.6f}'.format(*cell))

    print(f'\nNumber of atoms: {num_atoms}')
    print('\nAtomic positions (cart. coord.):' '\n x\t y\t z\t number of WF:')
    for pos, wan in zip(positions, wanniers):
        print('{:.6f} {:.6f} {:.6f} {:.6f}'.format(*pos, wan))

    print(
        f'\nTotal number of Wannier functions: {num_wannier_funcs}'
        f'\nPrint the complex part of the hamiltonian: {bool(print_complex)}')
    with open(args.output_file_name, 'w') as fp:
        print('Crystal structure of system:', file=fp)
        print('crystal axes: (cart. coord. )', file=fp)
        for cell in cell_vectors:
            print('{:.6f} {:.6f} {:.6f}'.format(*cell), file=fp)
        print('\nNumber of atoms: ', num_atoms, file=fp)
        print(
            '\nAtomic positions (cart. coord.):'
            '\n x\t y\t z\t number of WF:',
            file=fp)
        for pos, wan in zip(positions, wanniers):
            print('{:.6f} {:.6f} {:.6f} {:.6f}'.format(*pos, wan), file=fp)
        print('\nTotal number of Wannier functions: ',
              num_wannier_funcs,
              file=fp)
        print('\nPrint the complex part of the hamiltonian: ',
              print_complex,
              file=fp)
        
    # =======================================================================

    # move_x and move_y control the shift to the necessary block in Ham_R
    move_x = move_y = 0

    num_points = num_atoms * n_size.prod()
    index = np.zeros((num_points, 4), dtype=int)  # index = [i, j, k, atom]
    radius = np.zeros(num_points, dtype=float)

    # 'atom' is the index of the central atom
    with open(args.output_file_name, 'a') as fp:
        for atom, (pos, wan) in enumerate(zip(positions, wanniers)):
            # This function sorts atoms depending on the radius from the central atom with index 'atom'
            radius, index = coordination_sort(atom, num_atoms, n_min, n_max, cell_vectors,
                                              positions)

            print('-' * 69, file=fp)

            neighbor_num = 1
            sphere_num = 0
            for z, (rad, idx) in enumerate(zip(radius, index)):
                if z:
                    if abs(radius[z - 1] - rad) < 1e-4:
                        neighbor_num += 1
                    else:
                        neighbor_num = 1
                        sphere_num += 1

                if sphere_num == max_sphere_num:
                    break

                d3, at = idx[:3], idx[3]
                r = d3 @ cell_vectors + (positions[at] - pos)

                print(f'Hopping <a|H|b> between atom  {atom} (000)<-->atom  {at}',
                      '(',
                      *d3,
                      f') in sphere # {sphere_num}  '
                      'with radius  {:.6f}'.format(rad),
                      f' -- {neighbor_num}:',
                      file=fp)
                print('Radius vector is:  {:.6f} {:.6f} {:.6f}'.format(*r),
                      file=fp)

                move_y = wanniers[:at].sum()
                s = Ham_R[tuple(d3 + n_max)][move_x:,
                    move_y:][:wan, :wanniers[at]]

                fmt = '({0.real:.4f}  {0.imag:.4f}i)' if print_complex else '{0.real:.6f}'
                print('\n'.join('  '.join(fmt.format(item) for item in row)
                                for row in s),
                      file=fp)
                print(file=fp)

            move_x += wan

    print(f'This run was terminated on: {datetime.now()}')
    print(f'JOB DONE')
    print('=' * 69)



if __name__ == '__main__':
    main()

    
