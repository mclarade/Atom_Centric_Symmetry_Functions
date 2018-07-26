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


def file_list_builder():
    '''
    This function takes in the names of all of the files to be used and stores them in lists
    '''
    InputList = []
    for file in os.listdir(os.curdir):
        if file.endswith(args.InputExtension):
            InputList.append(file)
    InputList.sort()
    return InputList


def calculate_fcRij_matrix(distance_matrix):
    '''
    Calculates fcRij for the atom in question
    '''
    fcRij_matrix = 0.5 * (np.cos(np.pi * distance_matrix / args.Cutoff) + 1)
    fcRij_matrix[distance_matrix > args.Cutoff] = 0
    fcRij_matrix[distance_matrix == 0] = 0
    return fcRij_matrix


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
    G2 += math.exp(args.gausswidth * (distance - args.radial_distance) ** 2) * fcRij
    return G2


def calculate_g3(fcRij, distance):
    '''
    Takes in fcRij and distances, calculates and returns G3 for an atom, see (arxiv link)
    '''
    G3 = 0
    G3 += math.cos(distance * args.period_length) * fcRij
    return G3


def calculate_g4(fcRij, fcRik, fcRjk, distance_ab, distance_ac, distance_bc, angle):
    '''
    Calculates and returns G4 for an atom, takes in distances and list of angles see (arxiv link)
    '''
    G4 = 0
    G4 += ((1 + args.lambda_value * math.cos(angle)) ** args.angular_resolution
           * math.exp(-args.gausswidth * (distance_ab ** 2 + distance_ac ** 2 + distance_bc ** 2))
           * fcRij * fcRik * fcRjk) * 2 ** (1 - args.angular_resolution)
    return G4


def calculate_g5(fcRij, fcRik, distance_ab, distance_ac, angle):
#There are problems with the angle calculations
    '''
    Calculates and returns G5 for an atom, takes in distances and list of angles see (arxiv link)
    '''
    G5 = 0
    G5 += ((1 + args.lambda_value * math.cos(angle)) ** args.angular_resolution * math.exp(-args.gausswidth 
           * (distance_ab ** 2 + distance_ac ** 2)) * fcRij * fcRik
           * 2 ** (1 - args.angular_resolution))
    return G5


def calculate_angle(distance_ab, distance_ac, distance_bc):
    """
    called by detect_and_retrieve_angle to calculate angles
    """
    angle = math.acos((distance_ab ** 2 + distance_ac ** 2 - distance_bc ** 2) / (2 * distance_ab * distance_ac))
    return angle


def gen_energy_list(int_file):
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


def gen_symm_functions(matrix, labels, matrix_cutoff, array_dict_sorted_by_atom, atom_counter):
    """Calculates each individual symm function, and puts it into the appropriate cell of the numpy array"""
    symm_function_list = []
    for atom in labels:
        atom = filter(lambda x: not x.isdigit(), atom)
        if atom not in atom_counter.keys():
            atom_counter[atom] = 0
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
        atom_type = filter(lambda x: not x.isdigit(), labels[i])
        if atom_type not in array_dict_sorted_by_atom.keys():
            raise Exception('Unrecognized atom detected in database, please include it in the input dictionary using the -a or --atoms flag')
        for j in range(0, len(matrix)):
            if i == j:
                pass
            elif matrix_cutoff[i, j] < 1e-7:
                pass
            else:
                symm_function = []
                fcRij_matrix = calculate_fcRij_matrix(matrix)
                symm_function.append(matrix_cutoff[i, j])
                print i, j
                if args.G1flag == True:
                    G1_Total += calculate_g1(fcRij_matrix[i,j])
                if args.G2flag == True:
                    G2_Total += calculate_g2(fcRij_matrix[i,j], 
                                             matrix_cutoff[i, j])
                if args.G3flag == True:
                    G3_Total += calculate_g3(fcRij_matrix[i, j],
                                             matrix_cutoff[i, j])
                #start K loop here
                for k in range(0, len(matrix_cutoff)):
                    if i == k or j == k:
                       pass
                    elif matrix_cutoff[i, k] < 1e-7:
                        pass
                    else:
                        print i, j, k
                        if args.G4flag == True or args.G5flag == True:
                            angle = calculate_angle(
                                matrix_cutoff[i, j], matrix_cutoff[i, k], matrix[j, k])
                        if args.G4flag == True:
                            G4_Total += calculate_g4(fcRij_matrix[i, j], fcRij_matrix[i, k], fcRij_matrix[j, k], 
                                                        matrix_cutoff[i, j], matrix_cutoff[i, k], matrix[j, k], 
                                                        angle)
                        if args.G5flag == True:
                            G5_Total += calculate_g5(fcRij_matrix[i, j], fcRij_matrix[i, k], 
                                                        matrix_cutoff[i, j], matrix_cutoff[i, k], 
                                                        angle)
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
        print atom_counter[atom_type]
#        print array_dict_sorted_by_atom[atom_counter[atom_type], :]
        array_dict_sorted_by_atom[atom_type][atom_counter[atom_type],:] = symm_function_list
        print array_dict_sorted_by_atom[atom_type][atom_counter[atom_type], :]
        atom_counter[atom_type] += 1



def initialize_numpy_bins():
    '''Creates numpy bins for each atom type'''
    counter_dict = {}
    for atom_type in args.AtomInputList.keys():
        counter_dict.setdefault(atom_type, 0)
    OutputDimension2 = generate_output_dimensions()
    wavefunction_and_file_dict = {}
    InputList = file_list_builder()
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


def save_energy_data(energy_dict):
    for key in energy_dict.keys():
        energy_out = np.asarray(energy_dict[key])
        np.save(key+"_energy", energy_out)
        

def save_g_values(symm_dict):
    for key in symm_dict.keys():
        symm_out =  np.asarray(symm_dict[key])
        np.save(key+"_symm", symm_out)


def main(args):
    keylist, array_dict_sorted_by_atom = initialize_numpy_bins()
    energy_dict = {}
    symm_dict = {}
    for wavefunction in keylist:
        print wavefunction
        labels, distance_matrix, matrix_cutoff = retrieve_coordinates(wavefunction)
        symm_data = gen_symm_functions(distance_matrix, labels, matrix_cutoff, array_dict_sorted_by_atom, symm_dict)
        for intfile in wavefunction:
            atom_label, energy = gen_energy_list(intfile)
            if atom_label not in energy_dict:
                energy_dict[atom_label] = []
            energy_dict[atom_label].append(energy)
    save_energy_data(energy_dict)
    save_g_values(symm_dict)



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
    parser.add_argument('-c', '--cutoff',
                        dest='Cutoff',
                        help='Set cutoff distance in Bohr, default = 2.5, suggested values are between 2.0 and 11.0 Bohr',
                        type=float,
                        default=3)
    parser.add_argument('--G1',
                        dest='G1flag',
                        help='Set flag to calculate G1, default=True',
                        type=bool,
                        default=True)
    parser.add_argument('--G2',
                        dest='G2flag',
                        help='Set flag to calculate G2, default=False',
                        type=bool,
                        default=True)
    parser.add_argument('-g', '--gausswidth',
                        dest='gausswidth',
                        help='Set width of gaussian, in Bohr^-2, required for G2, G4 and G5, suggested values are between 0.01 and 5.00Bohr^-2',
                        type=float,
                        default=1.1)
    parser.add_argument('-d', '--radial_distance',
                        dest='radial_distance',
                        help='Set radial distance, in Bohr, required for G2, suggested values are between 2 and 10 Bohr',
                        type=float,
                        default=2.0)
    parser.add_argument('--G3',
                        dest='G3flag',
                        help='Set flag to calculate G3, default=False',
                        type=bool,
                        default=True)
    parser.add_argument('-p', '--period_length',
                        dest='period_length',
                        help='Set period length, in Bohr^-1, required for G3, suggested values are between 0.5 and 2.0 Bohr^-2',
                        type=float,
                        default=1.15)
    parser.add_argument('--G4',
                        dest='G4flag',
                        help='Set flag to calculate G4, default=False',
                        type=bool,
                        default=True)
    parser.add_argument('-l', '--lambda',
                        dest='lambda_value',
                        help='Set lambda, the maximum angle to 0 or pi radians by setting to +1 or -1, respectively, required for G4 and G5 default=1',
                        type=int,
                        default=-1)
    parser.add_argument('--G5',
                        dest='G5flag',
                        help='Set flag to calculate G5, default=False',
                        type=bool,
                        default=True)
    parser.add_argument('-r', '-resolution',
                        dest='angular_resolution',
                        help='Set angular resolution, required for G4 and G5, required to be int, suggested values are 1, 2, 4, 16 and 64 default=4',
                        type=int,
                        default=4)
    parser.add_argument('-a', '--atoms',
                        dest='AtomInputList',
                        help='Add to list of atoms to be inspected, takes input in the form Symbol:Name (eg, H:Hydrogen)',
                        type=dict,
                        default={'H':'Hydrogen','C':'Carbon','N':'Nitgrogen','O':'Oxygen'})
                        #add file reading logic
    args = parser.parse_args() 
    main(args)
#add in logic to throw error if G-functions are missing their required inputs
