#include <iostream>
#include <cmath>
#include <vector>
#include <array>

using namespace std;

array<double, 3> cross_product(array<double, 3> a, array<double, 3> b) {
    return {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]};
}

double norm(array<double, 3> a) {
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

array<double, 3> normalize(array<double, 3> a) {
    double n = norm(a);
    return {a[0]/n, a[1]/n, a[2]/n};
}

array<double, 3> rotate_point(array<double, 3> point, double angle, array<double, 3> line_direction_cosines) {
    angle = angle * M_PI / 180.0;
    array<array<double, 3>, 3> rotation_matrix = {{
        {cos(angle) + pow(line_direction_cosines[0], 2) * (1 - cos(angle)), line_direction_cosines[0] * line_direction_cosines[1] * (1 - cos(angle)) - line_direction_cosines[2] * sin(angle), line_direction_cosines[0] * line_direction_cosines[2] * (1 - cos(angle)) + line_direction_cosines[1] * sin(angle)},
        {line_direction_cosines[0] * line_direction_cosines[1] * (1 - cos(angle)) + line_direction_cosines[2] * sin(angle), cos(angle) + pow(line_direction_cosines[1], 2) * (1 - cos(angle)), line_direction_cosines[1] * line_direction_cosines[2] * (1 - cos(angle)) - line_direction_cosines[0] * sin(angle)},
        {line_direction_cosines[0] * line_direction_cosines[2] * (1 - cos(angle)) - line_direction_cosines[1] * sin(angle), line_direction_cosines[1] * line_direction_cosines[2] * (1 - cos(angle)) + line_direction_cosines[0] * sin(angle), cos(angle) + pow(line_direction_cosines[2], 2) * (1 - cos(angle))}
    }};

    return {rotation_matrix[0][0]*point[0] + rotation_matrix[0][1]*point[1] + rotation_matrix[0][2]*point[2], rotation_matrix[1][0]*point[0] + rotation_matrix[1][1]*point[1] + rotation_matrix[1][2]*point[2], rotation_matrix[2][0]*point[0] + rotation_matrix[2][1]*point[1] + rotation_matrix[2][2]*point[2]};
}

array<double, 3> get_fourth_atom_coordinates(array<double, 3> first_atom_coordinates, array<double, 3> second_atom_coordinates, array<double, 3> third_atom_coordinates, double bond_length, double bond_angle, double dihedral_angle) {
    array<double, 3> p = {second_atom_coordinates[0] - first_atom_coordinates[0], second_atom_coordinates[1] - first_atom_coordinates[1], second_atom_coordinates[2] - first_atom_coordinates[2]};
    array<double, 3> q = {third_atom_coordinates[0] - second_atom_coordinates[0], third_atom_coordinates[1] - second_atom_coordinates[1], third_atom_coordinates[2] - second_atom_coordinates[2]};

    array<double, 3> ABC_plane_normal = cross_product(p, q);
    ABC_plane_normal = normalize(ABC_plane_normal);

    array<double, 3> unit_vector_along_q = normalize(q);

    array<double, 3> initial_pos = {unit_vector_along_q[0] * bond_length, unit_vector_along_q[1] * bond_length, unit_vector_along_q[2] * bond_length};

    array<double, 3> rotated_pos = rotate_point(initial_pos, 180-bond_angle, ABC_plane_normal);

    array<double, 3> final_pos = rotate_point(rotated_pos, dihedral_angle, unit_vector_along_q);

    return {final_pos[0] + third_atom_coordinates[0], final_pos[1] + third_atom_coordinates[1], final_pos[2] + third_atom_coordinates[2]};
}

int main() {
    vector<array <double, 3> > atom_coordinates = {{{0, 0, 0}}, {{1.453, 0, 0}}, {{1.959, -0.628, -1.30}}};
    double bond_length = 1.3;
    double bond_angle = 120;
    double dihedral_angle = 30;

    array<double, 3> fourth_atom_coordinates = get_fourth_atom_coordinates(atom_coordinates[0], atom_coordinates[1], atom_coordinates[2], bond_length, bond_angle, dihedral_angle);

    cout << "Coordinates of fourth atom are: {" << fourth_atom_coordinates[0] << ", " << fourth_atom_coordinates[1] << ", " << fourth_atom_coordinates[2] << "}" << endl;

    return 0;
}
