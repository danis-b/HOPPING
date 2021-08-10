#include <math.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <nlohmann/json.hpp> //we use json format for input file

using json = nlohmann::json;


void
coordination_sort(int atom, int num_atoms, int n_min[3], int n_max[3], double cell_vectors[3][3], double **positions,
                  double *radius, int **index) {


    int i, j, k, x, z;
    int num_points;

    int n_size[3];
    double r[3];


    for (i = 0; i < 3; i++) {
        n_size[i] = n_max[i] - n_min[i] + 1; // Plus 1 for 0th
    }

    num_points = num_atoms * (n_size[0] * n_size[1] * n_size[2]);


    // Index filling
    i = n_min[0];
    j = n_min[1];
    k = n_min[2];
    x = 0;
    for (z = 0; z < num_points; z++) {
        index[z][0] = i;
        index[z][1] = j;
        index[z][2] = k;
        index[z][3] = x;
        r[0] = i * cell_vectors[0][0] + j * cell_vectors[1][0] + k * cell_vectors[2][0] +
               (positions[x][0] - positions[atom][0]);
        r[1] = i * cell_vectors[0][1] + j * cell_vectors[1][1] + k * cell_vectors[2][1] +
               (positions[x][1] - positions[atom][1]);
        r[2] = i * cell_vectors[0][2] + j * cell_vectors[1][2] + k * cell_vectors[2][2] +
               (positions[x][2] - positions[atom][2]);
        radius[z] = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        x++;
        if (x == num_atoms) {
            x = 0;
            k++;
            if (k == n_max[2] + 1) {
                k = n_min[2];
                j++;
            }
            if (j == n_max[1] + 1) {
                j = n_min[1];
                i++;
            }

        }
    }


    int Temp_Index[4];
    double Temp_Radius;

    //Simple bubble sort

    for (z = 0; z < num_points; z++) {
        for (x = z; x < num_points; x++) {
            if (radius[z] > radius[x]) {
                Temp_Index[0] = index[z][0];
                Temp_Index[1] = index[z][1];
                Temp_Index[2] = index[z][2];
                Temp_Index[3] = index[z][3];
                Temp_Radius = radius[z];

                index[z][0] = index[x][0];
                index[z][1] = index[x][1];
                index[z][2] = index[x][2];
                index[z][3] = index[x][3];
                radius[z] = radius[x];

                index[x][0] = Temp_Index[0];
                index[x][1] = Temp_Index[1];
                index[x][2] = Temp_Index[2];
                index[x][3] = Temp_Index[3];
                radius[x] = Temp_Radius;
            }
        }
    }


}


int main() {

    time_t td;
    int i, j, k, x, y, z;
    int n_size[3], n_min[3], n_max[3], vecs[3], orbs[2];
    int num_atoms, num_wannier_funcs, num_points, max_sphere_num, neighbor_num, sphere_num, move_x, move_y, print_complex, atom;

    int *wanniers, **index;

    double cell_vectors[3][3], ham_values[2], r[3];
    double **positions, *radius;
    std::complex<double> *****Ham_R;

    std::string input_file_name;


    td = time(NULL);

    std::cout << "Program HOPPING.x v.3.0 (c++) starts on ";
    std::cout << ctime(&td) << std::endl;
    std::cout << "=====================================================================" << std::endl;

    std::cout << "Enter the name of input file" << std::endl;
    std::cin >> input_file_name;


    std::ifstream hr;
    hr.open(input_file_name + "_hr.dat");
    if (!hr) {
        std::cout << "ERROR! Cannot open hopping file " << input_file_name << "_hr.dat" << std::endl;
        return 0;
    }



    //define the possible ranges for n_min and n_max from hr.dat file
    for (i = 0; i < 3; i++) {
        n_max[i] = n_min[i] = 0;
    }

    while (!hr.eof()) {
        hr >> vecs[0] >> vecs[1] >> vecs[2] >> orbs[0] >> orbs[1] >> ham_values[0] >> ham_values[1];

        for (i = 0; i < 3; i++) {
            if (n_max[i] < vecs[i]) n_max[i] = vecs[i];
            else if (n_min[i] > vecs[i]) n_min[i] = vecs[i];
        }

    }


    for (i = 0; i < 3; i++) {
        n_size[i] = n_max[i] - n_min[i] + 1; // Plus 1 for 0th
    }


    //===============================================================================
    //read information from input file in json format

    std::ifstream in;
    in.open(input_file_name + ".json");
    if (!in) {
        std::cout << "ERROR! Cannot open input file " << input_file_name << ".json" << std::endl;
        return 0;
    }


    json json_file;
    in >> json_file;

    in.close();


    // read int numbers from json file
    json_file.at("number_of_atoms").get_to(num_atoms);
    json_file.at("max_sphere_num").get_to(max_sphere_num);
    json_file.at("print_complex").get_to(print_complex);


    // create matrices for atomic position and number of wannier functions
    positions = new double *[num_atoms];
    for (i = 0; i < num_atoms; i++) {
        positions[i] = new double[3];
    }

    wanniers = new int[num_atoms];


    // read matrices from json file

    i = 0;
    for (auto &element : json_file["wanniers"]) {
        wanniers[i] = element;
        i++;
    }


    i = 0;
    for (auto &element : json_file["cell_vectors"]) {
        j = 0;
        for (auto &element1 : element) {
            cell_vectors[i][j] = element1;
            j++;
        }
        i++;
    }


    i = 0;
    for (auto &element : json_file["positions"]) {
        j = 0;
        for (auto &element1 : element) {
            positions[i][j] = element1;
            j++;
        }
        i++;
    }


    //total number of Wannier functions
    num_wannier_funcs = 0;
    for (i = 0; i < num_atoms; i++) {
        num_wannier_funcs += wanniers[i];
    }



    // Total hamiltonian
    Ham_R = new std::complex<double> ****[n_size[0]];
    for (i = 0; i < n_size[0]; i++) {
        Ham_R[i] = new std::complex<double> ***[n_size[1]];
        for (j = 0; j < n_size[1]; j++) {
            Ham_R[i][j] = new std::complex<double> **[n_size[2]];
            for (k = 0; k < n_size[2]; k++) {
                Ham_R[i][j][k] = new std::complex<double> *[num_wannier_funcs];
                for (x = 0; x < num_wannier_funcs; x++) {
                    Ham_R[i][j][k][x] = new std::complex<double>[num_wannier_funcs];
                }
            }
        }
    }


    hr.clear();
    hr.seekg(0, std::ios::beg);


    while (!hr.eof()) {
        hr >> vecs[0] >> vecs[1] >> vecs[2] >> orbs[0] >> orbs[1] >> ham_values[0] >> ham_values[1];

        Ham_R[vecs[0] + n_max[0]][vecs[1] + n_max[1]][vecs[2] + n_max[2]][orbs[0] - 1][orbs[1] -
                                                                                       1] = std::complex<double>(
                ham_values[0], ham_values[1]);
    }
    std::cout << "Hopping file " << input_file_name << "_hr.dat was  scanned  successfully" << std::endl;

    hr.close();



    //===============================================================================
    // print data to stdout and out-file
    

    std::ofstream out;
    out.open("HOPPING.dat");


    std::cout << "Crystal structure of system:" << std::endl;
    out << "Crystal structure of system:" << std::endl;

    std::cout << "crystal axes: (cart. coord. )" << std::endl;
    out << "crystal axes: (cart. coord. )" << std::endl;


    for (i = 0; i < 3; i++) {

        std::cout << std::fixed << std::setprecision(6) << cell_vectors[i][0] << "  " << cell_vectors[i][1] << "  "
                  << cell_vectors[i][2]
                  << std::endl;
        out << std::fixed << std::setprecision(6) << cell_vectors[i][0] << "  " << cell_vectors[i][1] << "  "
            << cell_vectors[i][2]
            << std::endl;
    }

    std::cout << std::endl;
    out << std::endl;

    std::cout << "Number of atoms: " << num_atoms << std::endl;
    out << "Number of atoms: " << num_atoms << std::endl;

    std::cout << std::endl << "Atomic positions (cart. coord.):\n x\t y\t z\t number of WF:" << std::endl;
    out << std::endl << "Atomic positions (cart. coord. ):\n x\t y\t z\t number of WF:" << std::endl;


    for (i = 0; i < num_atoms; i++) {
        std::cout << std::fixed << std::setprecision(6) << positions[i][0] << "  " << positions[i][1] << "  "
                  << positions[i][2] << "  "
                  << wanniers[i]
                  << std::endl;
        out << std::fixed << std::setprecision(6) << positions[i][0] << "  " << positions[i][1] << "  "
            << positions[i][2]
            << "  "
            << wanniers[i]
            << std::endl;
    }

    std::cout << "Total number of Wannier functions: " << num_wannier_funcs << std::endl;
    out << "Total number of Wannier functions: " << num_wannier_funcs << std::endl;


    std::cout << "Print the complex part of the hamiltonian: " << print_complex << std::endl;
    out << "Print the complex part of the hamiltonian: " << print_complex << std::endl;


    std::cout << std::endl;
    out << std::endl;





    //==============================================================================

    // move_x and move_y control the shift to the necessary block in Ham_R
    move_x = move_y = 0;

    num_points = num_atoms * (n_size[0] * n_size[1] * n_size[2]);


    // index = [i, j, k, atom]
    index = new int *[num_points];
    for (i = 0; i < num_points; i++) {
        index[i] = new int[4];
    }

    radius = new double[num_points];


    //'atom' is the index of the central atom
    for (atom = 0; atom < num_atoms; atom++) {

        //This function sorts atoms depending on the radius from the central atom with index 'atom'
        coordination_sort(atom, num_atoms, n_min, n_max, cell_vectors, positions, radius, index);


        out << "---------------------------------------------------------------------" << std::endl;

        neighbor_num = 1;
        sphere_num = 0;


        for (z = 0; z < num_points; z++) {

            if (z != 0) {
                if (fabs(radius[z - 1] - radius[z]) < 1E-4)neighbor_num++;
                else {
                    neighbor_num = 1;
                    sphere_num++;
                }
            }
            if (sphere_num == max_sphere_num) break;


            out << "Atom " << atom + 1 << "(000)<-->Atom " << index[z][3] + 1 << "(" << index[z][0] << index[z][1]
                << index[z][2] << ") in sphere #" << sphere_num << " with radius " << radius[z] << " is #"
                << neighbor_num
                << ":  " << std::endl;


            r[0] = index[z][0] * cell_vectors[0][0] + index[z][1] * cell_vectors[1][0] +
                   index[z][2] * cell_vectors[2][0] +
                   (positions[index[z][3]][0] - positions[atom][0]);
            r[1] = index[z][0] * cell_vectors[0][1] + index[z][1] * cell_vectors[1][1] +
                   index[z][2] * cell_vectors[2][1] +
                   (positions[index[z][3]][1] - positions[atom][1]);
            r[2] = index[z][0] * cell_vectors[0][2] + index[z][1] * cell_vectors[1][2] +
                   index[z][2] * cell_vectors[2][2] +
                   (positions[index[z][3]][2] - positions[atom][2]);

            out << std::fixed << std::setprecision(6) << "Radius vector is: " << r[0] << " " << r[1] << " " << r[2]
                << std::endl;

            move_y = 0;
            for (x = 0; x < index[z][3]; x++) {
                move_y += wanniers[x];
            }

            if (print_complex == 1) {
                for (x = move_x; x < wanniers[atom] + move_x; x++) {
                    for (y = move_y; y < wanniers[index[z][3]] + move_y; y++) {


                        out << std::fixed << std::setprecision(6)
                            << Ham_R[index[z][0] + n_max[0]][index[z][1] + n_max[1]][index[z][2] + n_max[2]][x][y]
                            << "  ";

                    }
                    out << std::endl;
                }
                out << std::endl;
            } else {
                for (x = move_x; x < wanniers[atom] + move_x; x++) {
                    for (y = move_y; y < wanniers[index[z][3]] + move_y; y++) {


                        out << std::fixed << std::setprecision(6)
                            << (Ham_R[index[z][0] + n_max[0]][index[z][1] + n_max[1]][index[z][2] +
                                                                                      n_max[2]][x][y]).real()
                            << "  ";

                    }
                    out << std::endl;
                }
                out << std::endl;
            }

        }
        move_x += wanniers[atom];
    }


    td = time(NULL);
    std::cout << std::endl << "This run was terminated on: " << std::endl;
    std::cout << ctime(&td) << std::endl;

    std::cout << "JOB DONE" << std::endl;
    std::cout << "=====================================================================" << std::endl;


    delete[] Ham_R;
    delete[] radius;
    delete[] index;
    delete[] positions;
    delete[] wanniers;


}