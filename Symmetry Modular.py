from __future__ import division
import os
import math
import numpy as np
import argparse
from itertools import combinations

from scipy.spatial.distance import pdist, squareform

#adapt this
CarbonWritten = 0
OxygenWritten = 0
HydrogenWritten = 0
HEnergyList = []
CEnergyList = []
OEnergyList = []


def list_builder():
    '''
    This function takes in the names of all of the files to be used and stores them in lists
    '''
    InputList = []
    WaveFunctionList = []
    for file in os.listdir(os.curdir):
        if file.endswith(args.InputExtension):
            InputList.append(file)
        if file.endswith(args.WaveFunctionExtension):
            WaveFunctionList.append(file)
    InputList.sort()
    WaveFunctionList.sort()
    return InputList, WaveFunctionList


def calculate_fcRij(distance):
    '''
    Calculates fcRij for the atom in question
    '''
    if distance <= 1e-7:
        pass
    if distance > args.Cutoff:
        fcRij = 0
    else:
        fcRij = 0.5 * (math.cos(math.pi * distance / args.Cutoff) + 1)
    return fcRij


# redundant, but left in for clarity
def calculate_g1(fcRij):
    '''
    Sums up fcRij to make G1 for an atom, see (arxiv link)
    '''
    G1 = 0
    G1 += fcRij
    return G1


def calculate_g2(fcRij, distance):
    '''
    Takes in fcRij and distances, calculates and returns G2 for an atom, see (arxiv link)
    '''
    G2 = 0
    args.gausswidth
    G2 += math.exp(gausswidth * (distance) ** 2) * fcRij
    return G2


def calculate_g3(fcRij, distance):
    '''
    Takes in fcRij and distances, calculates and returns G3 for an atom, see (arxiv link)
    '''
    G3 = 0
    global period_length
    G3 += math.cos(distance * period_length) * fcRij
    return G3


def calculate_g4(fcRij, distance, angle):
    '''
    Calculates and returns G4 for an atom, takes in distances and list of angles see (arxiv link)
    '''
    G4 = 0
    G4 += (1 + args.Lamba * math.cos(angle)) ** angular_resolution * math.exp(-gausswidth * (distance12 ** 2 + distance13 ** 2 + distance23 ** 2)) * calculate_fcRij(distance12) * calculate_fcRij(distance13) * calculate_fcRij(distance23)
    G4 = G4 * 2 ** (1 - angular_resolution)
    return G4


def calculate_g5(fcRij, distance, angle):
#There are problems with the angle calculations
    '''
    Calculates and returns G5 for an atom, takes in distances and list of angles see (arxiv link)
    '''
    G5 = 0
    global gausswidth, Lamba, angular_resolution
    G5 += (1 + Lamba * math.cos(angle)) ** angular_resolution * math.exp(-gausswidth * (distance12 ** 2 + distance13 ** 2)) * calculate_fcRij(distance12) * calculate_fcRij(distance13)
    G5 = G5 * 2 ** (1 - angular_resolution)
    return G5


def detect_and_retrieve_angle(atomic_coordinates_cutoff):
    """
    Takes in 3d coordinates, calculates and returns every angle between every atom that falls under the cutoff
    """
    neighbors_list = []
    angle_list = []
    for i, element in enumerate(atomic_coordinates_cutoff):
        if element >= 1e-7:
            neighbors_list.append(element)
    neighbor_combinations = combinations(neighbors_list, 2)
    for combination in neighbor_combinations:
        distance3 = atomic_coordinates_cutoff[combo]
        angle = calculate_angle(combination, distance3)
        angle_list.append(angle)
    return angle_list


def calculate_angle(dists_from_central_atom, distance_between_other_atoms):
    """
    called by detect_and_retrieve_angle to calculate angles
    """
    angle = math.acos((dists_from_central_atom[0] ** 2 + dists_from_central_atom[1] ** 2 - distance_between_other_atoms ** 2) / (2 * dists_from_central_atom[0] * dists_from_central_atom[1]))
    return angle


def gen_energy_list(int_file, array_dict):
    '''
    Reads energy values out of every .int file and stores them into a list to be saved
    '''
    with open(int_file, 'r') as atomic_file:
        atomic_lines = atomic_file.readlines()
        for line in atomic_lines:
            if line == atomic_lines[6]:
                atom_label = line[0]
            if line.startswith('              K'):
                floatenergy = float(line.split()[3])
        return atom_label, floatenergy


def retrieve_coordinates(wavefunction):
    '''
    Extracts coordinates from wavefunction file, and generates a distance matrix and cutoff distance matrix
    '''
    print wavefunction[0]
    root_of_filename = wavefunction[0].split('.')[0]
    filename = root_of_filename+'.wfn'
    with open(filename, 'r') as waveinput:
        wavelines = waveinput.readlines()
        label_list = []
        x_list = []
        y_list = []
        z_list = []
        for line in wavelines:
            linesplit = line.split()
            if (line == wavelines[0]) or (line == wavelines[1]):
                pass
            elif line.startswith("CENTRE ASSIGNMENTS"):
                break
            else:
                label_list.append((str(linesplit[0]) + str(linesplit[1])))
                x_list.append(float(linesplit[4]))
                y_list.append(float(linesplit[5]))
                z_list.append(float(linesplit[6]))
    position_matrix = np.stack((x_list, y_list, z_list), axis=-1)
    distance_list = pdist(position_matrix)
    distance_matrix = squareform(distance_list)
    matrix_cutoff = np.copy(distance_matrix)
    matrix_cutoff[matrix_cutoff >= args.Cutoff] = 0
    return label_list, distance_matrix, matrix_cutoff


def generate_output_dimensions():
    """Sets the size of second dimension for numpy output"""
    OutputDimension2 = 0 
    if args.G1flag:
        OutputDimension2 += 1
    if args.G2flag:
        OutputDimension2 += 1
    if args.G3flag:
        OutputDimension2 += 1
    if args.G4flag:
        OutputDimension2 += 1
    if args.G5flag:
        OutputDimension2 += 1
    return OutputDimension2


def gen_symm_functions(matrix, labels, matrix_cutoff):
    """Calculates each individual symm function, and puts it into the appropriate cell of the numpy array"""
    symm_function_list = []
    if args.G1flag == True:
        G1_Total = 0
    if args.G2flag == True:
        G2_Total = 0
    if args.G3flag == True:
        G3_Total = 0
    if args.G4flag == True:
        G4_Total = 0
    if args.G5flag == True:
        G5_Total = 0
    for i in range(0, len(matrix)):
        fcRij_list = []
        symm_function = []
        for distance in matrix[i]:
            if distance > 1e-7:
                symm_function.append(distance)
                fcRij_list.append(calculate_fcRij(distance))
            atom_label = labels[i][0]
        for fcRij in fcRij_list:
            if args.G1flag == 'Y':
                G1_Total += calculate_g1(fcRij)
            if args.G2flag == 'Y':
                G2_Total += calculate_g2(fcRij, distance)
            if args.G3flag == 'Y':
                G3_Total += calculate_g3(fcRij, distance)
            if args.G4flag == 'Y':
                angles = detect_and_retrieve_angle(matrix_cutoff)
                for angle in angles:
                    G4_Total += calculate_g4(fcRij, distance, angle)
            if args.G5flag == 'Y':
                angles = detect_and_retrieve_angle(matrix_cutoff)
                for angle in angles:
                    G5_Total += calculate_g5(fcRij, distance, angle) 
    if args.G1flag == True:
        symm_function_list.append(G1_Total)
    if args.G2flag == True:
        symm_function_list.append(G2_Total)
    if args.G3flag == True:
        symm_function_list.append(G3_Total)
    if args.G4flag == True:
        symm_function_list.append(G4_Total)
    if args.G5flag == True:
        symm_function_list.append(G5_Total)
    return symm_function_list

def initialize_numpy_bins():
    '''Creates numpy bins for each atom type'''
    counter_dict = {}
    for atom_type in args.AtomInputList.keys():
        counter_dict.setdefault(atom_type, 0)
    OutputDimension2 = generate_output_dimensions()
    wavefunction_and_file_dict = {}
    InputList, WaveFunctionList = list_builder()
    for atomfilename in InputList:
        wavefunction = atomfilename.split(args.WaveFunctionExtension)
        if wavefunction[0] not in wavefunction_and_file_dict:
            wavefunction_and_file_dict[wavefunction[0]] = [atomfilename]
        else:
            #This may be redundant
            atomfilelist = wavefunction_and_file_dict[wavefunction[0]]
            atomfilelist.append(atomfilename)
            wavefunction_and_file_dict[wavefunction[0]] = atomfilelist
        for atom_type in args.AtomInputList.keys():
            string_check = args.WaveFunctionExtension + atom_type
            if string_check in atomfilename:
                counter_dict[atom_type] += 1
                break
    array_dict = {}
    for atom_type in args.AtomInputList.keys():
        if counter_dict[atom_type] == 0:
            pass
        else:
            #look here for output dimensions
            dimension0 = counter_dict[atom_type]
            array_dict[atom_type] = np.zeros((dimension0, OutputDimension2))
    keylist = wavefunction_and_file_dict.values()
    keylist.sort()
    return keylist, array_dict


def save_data(energy_dict):
    for key in energy_dict.keys():
        energy_out = np.asarray(energy_dict[key])
        np.save(key+"_energy",energy_out)
        

def main(args):
    keylist, array_dict_sorted_by_atom = initialize_numpy_bins()
    energy_dict = {}
    for wavefunction in keylist:
        labels, distance_matrix, matrix_cutoff = retrieve_coordinates(wavefunction)
        gen_symm_functions(distance_matrix, labels, matrix_cutoff)
        for intfile in wavefunction:
            atom_label, energy = gen_energy_list(intfile, array_dict_sorted_by_atom)
            if atom_label not in energy_dict:
                energy_dict[atom_label] = []
            energy_dict[atom_label].append(energy) 
    save_data(energy_dict)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='This script automates the conversion of wfn and int files into atom-centric symmetry functions for use with neural network inputs. Note that each wavefunction must have a unique name, or behavior may become unpredictable')
    parser.add_argument('-x', '--extension',
                        dest='InputExtension',
                        help='Extension of aimpac output scripts',
                        default='.int')
    parser.add_argument('-w', '--wavefunction',
                        dest='WaveFunctionExtension',
                        help='Extension of wavefunction files',
                        default='.wfn')
    parser.add_argument('--G1',
                        dest='G1flag',
                        help='Set flag to calculate G1, default=True',
                        type=bool,
                        default=True)
    parser.add_argument('--G2',
                        dest='G2flag',
                        help='Set flag to calculate G2, default=False',
                        type=bool,
                        default=False)
    parser.add_argument('--G3',
                        dest='G3flag',
                        help='Set flag to calculate G3, default=False',
                        type=bool,
                        default=False)
    parser.add_argument('--G4',
                        dest='G4flag',
                        help='Set flag to calculate G4, default=False',
                        type=bool,
                        default=False)
    parser.add_argument('--G5',
                        dest='G5flag',
                        help='Set flag to calculate G5, default=False',
                        type=bool,
                        default=False)
    parser.add_argument('-r', '-resolution',
                        dest='AngularResolution',
                        help='Set angular resolution, default=5.0',
                        type=float,
                        default=5.0)
    parser.add_argument('-g', '--gausswidth',
                        dest='GaussWidth',
                        help='Set width of gaussian',
                        type=float,
                        default=3.0)
    parser.add_argument('-l', '--lambda',
                        dest='lambda',
                        help='Set lambda, default=1',
                        type=int,
                        default=1)
    parser.add_argument('-c', '--cutoff',
                        dest='Cutoff',
                        help='Set cutoff distance in bohr, defualt = 4.5',
                        type=float,
                        default=4.5)
    parser.add_argument('-a', '--atoms',
                        dest='AtomInputList',
                        help='Add to list of atoms to be inspected, takes input in the form Symbol:Name (eg, Li:Lithium)',
                        type=dict,
                        default={'H':'Hydrogen','C':'Carbon','O':'Oxygen'})
                        #add file reading logic
    args = parser.parse_args() 
    main(args)
