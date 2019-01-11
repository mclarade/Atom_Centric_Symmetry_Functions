from __future__ import division
import os
import math
import numpy as np
import argparse
import json
from itertools import combinations

from scipy.spatial.distance import pdist, squareform


def str2bool(string):
    if string.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif string.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

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


def calculate_fcRij_matrix(distance_matrix, cutoff):
    '''
    Calculates fcRij for the atom in question
    '''
    fcRij_matrix = 0.5 * (np.cos(np.pi * distance_matrix / cutoff) + 1)
    fcRij_matrix[distance_matrix > cutoff] = 0
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


def calculate_g2(fcRij, distance, params):
    '''
    Takes in fcRij and distances, calculates and returns G2 for an atom, see (arxiv link)
    '''
    G2 = 0
    G2 += math.exp(params['gausswidth'] * (distance - params["radial_distance"]) ** 2) * fcRij
    return G2


def calculate_g3(fcRij, distance, params):
    '''
    Takes in fcRij and distances, calculates and returns G3 for an atom, see (arxiv link)
    '''
    G3 = 0
    G3 += math.cos(distance * params["period_length"]) * fcRij
    return G3


def calculate_g4(fcRij, fcRik, fcRjk, distance_ab, distance_ac, distance_bc, angle, params):
    '''
    Calculates and returns G4 for an atom, takes in distances and list of angles see (arxiv link)
    '''
    G4 = 0
    G4 += ((1 + params["lambda_value"] * math.cos(angle)) ** params["angular_resolution"]
           * math.exp(-params["gausswidth"] * (distance_ab ** 2 + distance_ac ** 2 + distance_bc ** 2))
           * fcRij * fcRik * fcRjk) * 2 ** (1 - params["angular_resolution"])
    return G4


def calculate_g5(fcRij, fcRik, distance_ab, distance_ac, angle, params):
#There are problems with the angle calculations
    '''
    Calculates and returns G5 for an atom, takes in distances and list of angles see (arxiv link)
    '''
    G5 = 0
    G5 += ((1 + params["lambda_value"] * math.cos(angle)) ** params["angular_resolution"] * math.exp(-params["gausswidth"] 
           * (distance_ab ** 2 + distance_ac ** 2)) * fcRij * fcRik
           * 2 ** (1 - params["angular_resolution"]))
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


def retrieve_coordinates(wavefunction, params):
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
#TODO: Make masks instead of matricies, and pass the masks around
    position_matrix = np.stack((x_list, y_list, z_list), axis=-1)
    distance_list = pdist(position_matrix)
    distance_matrix = squareform(distance_list)
    matrix_cutoff = np.copy(distance_matrix)
    cutoff_matricies = {}
    for cutoff_value in params["cutoff"]:
        cutoff_matricies[cutoff_value] = matrix_cutoff[matrix_cutoff >= cutoff_value] = 0
    return label_list, distance_matrix, cutoff_matricies


def generate_output_dimensions(params):
    """Sets the size of second dimension for numpy output"""
    OutputDimension2 = 0 
    if args.G1flag:
        OutputDimension2 += len(params["cutoff"])
    if args.G2flag:
        OutputDimension2 += len(params["cutoff"]) * len(params["gausswidth"]) * len(params["radial_distance"])
    if args.G3flag:
        OutputDimension2 += len(params["cutoff"]) * len(params["period_length"])
    if args.G4flag:
        OutputDimension2 += len(params["cutoff"]) * len(params["gausswidth"]) * len(params["lambda_value"]) * len(params["angular_resolution"])
    if args.G5flag:
        OutputDimension2 += len(params["cutoff"]) * len(params["gausswidth"]) * len(params["lambda_value"]) * len(params["angular_resolution"])
    return OutputDimension2


#TODO: Make this work with lists
def gen_symm_functions(matrix, labels, cutoff_matricies, array_dict_sorted_by_atom, atom_counter, params):
    """Calculates each individual symm function, and puts it into the appropriate cell of the numpy array"""
    for atom in labels:
        atom = filter(lambda x: not x.isdigit(), atom)
        if atom not in atom_counter.keys():
            atom_counter[atom] = 0
    for i in range(0, len(matrix)):
        if args.G1flag == True:
            G1_Total = []
        if args.G2flag == True:
            G2_Total = {}
        if args.G3flag == True:
            G3_Total = {}
        if args.G4flag == True:
            G4_Total = {}
        if args.G5flag == True:
            G5_Total = {}
        symm_function_list = []
        atom_type = filter(lambda x: not x.isdigit(), labels[i])
        if atom_type not in array_dict_sorted_by_atom.keys():
            raise Exception('Unrecognized atom detected in database, please include it in the input dictionary using the -a or --atoms flag')
        for j in range(0, len(matrix)):
            for cutoff_value in params["cutoff"]:
                if i == j:
                    pass
                elif cutoff_matricies[cutoff_value][i, j] < 1e-7:
                    pass
                else:
                    fcRij_matrix = calculate_fcRij_matrix(matrix, cutoff_value)
                    if args.G1flag == True:
                        for cutoff_value in params["cutoff"]:
                            G1_Total += calculate_g1(
                                fcRij_matrix[i, j])
                    if args.G2flag == True:
                        G2_Total += calculate_g2(fcRij_matrix[i,j], 
                                                cutoff_matricies[cutoff_value][i, j], params )
                    if args.G3flag == True:
                        G3_Total += calculate_g3(fcRij_matrix[i, j],
                                                cutoff_matricies[cutoff_value][i, j], params)
                    #start K loop here
                    for k in range(0, len(cutoff_matricies[cutoff_value])):
                        if i == k or j == k:
                            pass
                        elif cutoff_matricies[cutoff_value][i, k] < 1e-7:
                            pass
                        else:
                            if args.G4flag == True or args.G5flag == True:
                                angle = calculate_angle(
                                    cutoff_matricies[cutoff_value][i, j], cutoff_matricies[cutoff_value][i, k], matrix[j, k])
                            if args.G4flag == True:
                                G4_Total += calculate_g4(fcRij_matrix[i, j], fcRij_matrix[i, k], fcRij_matrix[j, k], 
                                                            cutoff_matricies[cutoff_value][i, j], cutoff_matricies[cutoff_value][i, k], matrix[j, k], 
                                                         angle, params)
                            if args.G5flag == True:
                                G5_Total += calculate_g5(fcRij_matrix[i, j], fcRij_matrix[i, k], 
                                                            cutoff_matricies[cutoff_value][i, j], cutoff_matricies[cutoff_value][i, k], 
                                                         angle, params)
        array_dict_sorted_by_atom[atom_type][atom_counter[atom_type],:] = symm_function_list
        atom_counter[atom_type] += 1
    return array_dict_sorted_by_atom


def initialize_numpy_bins(AtomInputList, OutputDimension2):
    '''Creates numpy bins for each atom type'''
    counter_dict = {}
    for atom_type in AtomInputList.keys():
        counter_dict.setdefault(atom_type, 0)
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
        for atom_type in AtomInputList.keys():
            string_check = args.WaveFunctionExtension + atom_type
            if string_check in atomfilename:
                counter_dict[atom_type] += 1
                break
    array_dict = {}
    for atom_type in AtomInputList.keys():
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
        

def save_g_values(symm_data):
    for key in symm_data.keys():
        symm_out =  np.asarray(symm_data[key])
        np.save(key+"_symm", symm_out)


def main(args):
    with open('Atom_Dict.json') as AtomsIn:
        AtomInputList = json.loads(AtomsIn.read())
    with open('Parameters.json') as params_in:
        params = json.loads(params_in.read())
    OutputDimension2 = generate_output_dimensions(params)
    keylist, array_dict_sorted_by_atom = initialize_numpy_bins(AtomInputList, OutputDimension2)
    energy_dict = {}
    counter_dict = {}
    for wavefunction in keylist:
        print wavefunction
        labels, distance_matrix, cutoff_matricies = retrieve_coordinates(wavefunction, params)
        array_dict_sorted_by_atom = gen_symm_functions(
            distance_matrix, labels, cutoff_matricies, array_dict_sorted_by_atom, counter_dict, params)
        for intfile in wavefunction:
            atom_label, energy = gen_energy_list(intfile)
            if atom_label not in energy_dict:
                energy_dict[atom_label] = []
            energy_dict[atom_label].append(energy)
    save_energy_data(energy_dict)
    save_g_values(array_dict_sorted_by_atom)


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
                        type=str2bool,
                        default=True)
    parser.add_argument('--G2',
                        dest='G2flag',
                        help='Set flag to calculate G2, default=False',
                        type=str2bool,
                        default=True)
    parser.add_argument('--G3',
                        dest='G3flag',
                        help='Set flag to calculate G3, default=False',
                        type=str2bool,
                        default=True)
    parser.add_argument('--G4',
                        dest='G4flag',
                        help='Set flag to calculate G4, default=False',
                        type=str2bool,
                        default=True)
    parser.add_argument('--G5',
                        dest='G5flag',
                        help='Set flag to calculate G5, default=False',
                        type=str2bool,
                        default=True)
    args = parser.parse_args() 
    main(args)
#add in logic to throw error if G-functions are missing their required inputs
