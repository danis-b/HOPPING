# HOPPING.x

This script converts seedname_hr.dat hamiltonian from  the wannier90 package to block-divided form between orbitals and prints these blocks depending on the geometry given by seedname.json file. Python version script can be used as is, while c++ version needs to be compiled with the additional [json](https://github.com/nlohmann/json) library.

# Usage 

As an example let's consider PBE band structure of Cu2GeO4 system [[Phys. Rev. B 100, 214401 (2019)](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.100.214401)], which was parametrized by Wannier functions based on Cu(d), O(p) and Ge(s) initial projectors. Total number of Wannier functions is 46, corresponding to 20(5x4) *d* orbitals of copper, 24(3x8) *p* orbitals of oxygen and 2(1x2) *s* orbitals of germanium atoms:

![alt text](https://github.com/danis-b/HOPPING/blob/main/example/bands.png)

seedname.json file contains the following information:

```json
   "cell_vectors": [[5.596748, 0.000000,  0.000000], [0.000000,  5.596748,  0.000000], [2.798374, 2.798374, 4.700648]],
   "number_of_atoms": 14,
   "max_sphere_num": 10,
   "print_complex": 0,
   "positions": [[2.79837,   4.19756,   1.17516], [6.99593,   5.59675,   3.52549], [2.79837,   1.39919,   1.17516], [4.19756,   5.59675,   3.52549], [5.59675,   6.94556,   3.37506], [1.44956,   5.59675,   1.02474], [5.59675,   4.24793,   3.37506], [4.14719,   5.59675,   1.02474], [2.79837,   4.24793,   3.67591], [1.44956,   2.79837,   1.32558], [2.79837,   6.94556,   3.67591], [4.14719,   2.79837,   1.32558], [0.00000,   0.00000,   0.00000], [5.59675,   2.79837,   2.35032]],
   "wanniers": [5,5,5,5,3,3,3,3,3,3,3,3,1,1]
```
* cell_vectors - unit cell vectors;
* number_of_atoms - total number of atoms, parametrized by Wannier functions;
* positions - cartesian coordinates of these atoms;
* wanniers - number of Wannier functions of these atoms;
* max_sphere_num - maximum number of coordination sphere to print the hamiltonian blocks between pair of atoms;
* print_complex - print the imaginary part of hamiltonian blocks. 

Both version of script needs to be started at the same folder with seedname.json and seedname_hr.dat file. **Please, make sure that the additional lines before hopping parameters in seedname_hr.dat file are removed.**




