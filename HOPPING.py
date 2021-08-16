import numpy as np
import json
from datetime import datetime


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


if __name__ == '__main__':
    print("Program HOPPING.x v.3.0 (python) starts on ", datetime.now())
    print('=' * 69)
    input_file_name = input('Enter the name of input file (name.json and name_hr.dat)\n')

    with open(f'{input_file_name}_hr.dat') as fp:
        rows = (line.split() for line in fp)
        data = [([int(u) for u in row[:3]], [int(u) for u in row[3:5]],
                 [float(u) for u in row[5:]]) for row in rows]

    # [N, 3] vectors, [N, 2] orbitals, [N, 2] hamiltonian
    vecs, orbs, ham_values = map(np.array, zip(*data))
    ham_values = ham_values.astype('f8').view(
        'c16').ravel()  # View as complex [N]

    n_min, n_max = vecs.min(0), vecs.max(0)  # [3]
    n_size = n_max - n_min + 1  # Plus 1 for 0th

    # =======================================================================
    # Read information from input file in json format
    with open(f'{input_file_name}.json') as fp:
        data = json.load(fp)

    num_atoms = data['number_of_atoms']  # M
    positions = np.array(data['positions'])  # array[M, 3]
    wanniers = np.array(data['wanniers'])  # array[M]
    cell_vectors = np.array(data['cell_vectors'])  # array[3, 3]
    max_sphere_num = data['max_sphere_num']
    print_complex = data['print_complex']  # Pretty repr

    # Some checks
    assert num_atoms == len(wanniers) == len(positions)
    num_wannier_funcs = wanniers.sum()

    # Total hamiltonian
    Ham_R = np.zeros((*n_size, num_wannier_funcs, num_wannier_funcs),
                   dtype='c16')
    Ham_R[(*(vecs + n_max).T, *(orbs.T - 1))] = ham_values

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
    with open('HOPPING.dat', 'w') as fp:
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
    with open('HOPPING.dat', 'a') as fp:
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

                print(f'Atom  {atom + 1} (000)<-->Atom  {at + 1}',
                      '(',
                      *d3,
                      f') in sphere # {sphere_num}  '
                      'with radius  {:.6f}'.format(rad),
                      f' is # {neighbor_num} :',
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
