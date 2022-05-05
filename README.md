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

As a result, one can find HOPPING.dat file, which contains the information from seedname.json file and output data:

==================================================================

Atom 0(000)<-->Atom 0(000) in sphere #0 with radius 0.000000 is #1:  
Radius vector is: 0.000000 0.000000 0.000000 \
7.184094  0.000004  0.051431  -0.019958  0.000000  
0.000004  6.783953  -0.000000  -0.000013  -0.050932  
0.051431  -0.000000  6.735947  -0.038462  0.000001  
-0.019958  -0.000013  -0.038462  6.223730  0.000002  
0.000000  -0.050932  0.000001  0.000002  7.065973  

Atom 0(000)<-->Atom 5(000) in sphere #1 with radius 1.949268 is #1:  
Radius vector is: -1.348810 1.399190 -0.150420 \
0.137210  0.318561  -0.327827  
0.454855  -0.128865  0.153048  
-0.481865  0.123510  -0.066666  
0.019533  0.491217  0.435620  
0.234957  0.849009  -0.890015  

........
==================================================================

For example, the last table gives the hamiltonian block 5x3 between Atom 0 (Cu(*d*)) and Atom 5 (O(*p*)) with the distance between them 1.949268 Å and radius vector **R** = (-1.348810, 1.399190, -0.150420)Å;




