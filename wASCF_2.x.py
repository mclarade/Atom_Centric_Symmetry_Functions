from __future__ import division
import os
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


def label_to_charge(label, AtomInputList):
    """
    label_to_charge takes in atom label and returns atom charge
    """
    if label in AtomInputList.keys():
        return AtomInputList[label][1]


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
    if fcRij_matrix.shape == ():
        return fcRij_matrix
    fcRij_matrix[distance_matrix > cutoff] = 0
    fcRij_matrix[distance_matrix == 0] = 0
    return fcRij_matrix


def radial_weight_calculator(fcRij, distance, atom_charge_j, gausswidth):
    rad_weight = (atom_charge_j) * np.exp(-gausswidth * (distance - gausswidth) ** 2) * fcRij
    return rad_weight


def angular_weight_calculator(fcRij, fcRik, fcRjk, distance_ij, distance_ik, distance_jk, atom_charge_j, atom_charge_k, gausswidth, angular_resolution, lambda_value, alternative_charge_calc = false):
    #TODO: think of better vairable name than charge_calc
    if alternative_charge_calc:
        charge_calc = atom_charge_j * atom_charge_k / (atom_charge_j + atom_charge_k)
    else:
        charge_calc = atom_charge_j * atom_charge_k
    ang_weight = ((1 + lambda_value * np.cos(angle)) ** angular_resolution * 
                charge_calc * 2 ** (1 - angular_resolution) * 
                np.exp(-gausswidth*(distance_ij-angular_resolution)) * 
                np.exp(-gausswidth*(distance_ik-angular_resolution)) *
                np.exp(-gausswidth*(distance_jk-angular_resolution)) *
                fcRij * fcRik * fcRjk)

# redundant, but left in for clarity
def calculate_g1(fcRij):
    '''
    Sums up fcRij to make G1 for an atom, see (arxiv link)
    '''
    G1 = fcRij
    return G1


def calculate_g2(fcRij, distance, gausswidth, radial_distance):
    '''
    Takes in fcRij and distances, calculates and returns G2 for an atom, see (arxiv link)
    '''
    G2 = np.exp(-gausswidth * (distance - radial_distance) ** 2) * fcRij
    return G2


def calculate_g3(fcRij, distance, period_length):
    '''
    Takes in fcRij and distances, calculates and returns G3 for an atom, see (arxiv link)
    '''
    G3 = 0
    G3 += np.cos(distance * period_length) * fcRij
    return G3


def calculate_g4(fcRij, fcRik, fcRjk, distance_ij, distance_ik, distance_jk, angle, lambda_value, angular_resolution, gausswidth):
    '''
    Calculates and returns G4 for an atom, takes in distances and list of angles see (arxiv link)
    '''
    G4 = 0
    G4 += ((1 + lambda_value * np.cos(angle)) ** angular_resolution
           * np.exp(-gausswidth * (distance_ij ** 2 + distance_ik ** 2 + distance_jk ** 2))
           * fcRij * fcRik * fcRjk) * 2 ** (1 - angular_resolution)
    return G4


def calculate_g5(fcRij, fcRik, distance_ij, distance_ik, angle, lambda_value, angular_resolution, gausswidth):
#There are problems with the angle calculations
    '''
    Calculates and returns G5 for an atom, takes in distances and list of angles see (arxiv link)
    '''
    G5 = 0
    G5 += ((1 + lambda_value * np.cos(angle)) ** angular_resolution * 
            np.exp(-gausswidth * (distance_ij ** 2 + distance_ik ** 2)) * fcRij * fcRik
            * 2 ** (1 - angular_resolution))
    return G5


def calculate_angle(distance_ij, distance_ik, distance_jk):
    """
    called by detect_and_retrieve_angle to calculate angles
    """
    angle = np.arccos((distance_ij ** 2 + distance_ik ** 2 - distance_jk ** 2) / (2 * distance_ij * distance_ik))
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
    position_matrix = np.stack((x_list, y_list, z_list), axis=-1)
    distance_list = pdist(position_matrix)
    distance_matrix = squareform(distance_list)
    return label_list, distance_matrix


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


def initialize_numpy_bins(AtomInputList, OutputDimension2):
    '''Creates numpy bins for each atom type'''
    total_dict = {}
    counter_dict = {}
    for atom_type in AtomInputList.keys():
        total_dict.setdefault(atom_type, 0)
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
                total_dict[atom_type] += 1
                break
    array_dict = {}
    for atom_type in AtomInputList.keys():
        if total_dict[atom_type] == 0:
            pass
        else:
            #look here for output dimensions
            dimension0 = total_dict[atom_type]
            array_dict[atom_type] = np.zeros((dimension0, OutputDimension2))
    keylist = wavefunction_and_file_dict.values()
    keylist.sort()
    return keylist, array_dict, counter_dict


def save_energy_data(energy_dict):
    for key in energy_dict.keys():
        energy_out = np.asarray(energy_dict[key])
        np.save(key+"_energy", energy_out)
        print key, energy_out

def save_g_values(symm_data):
    for key in symm_data.keys():
        symm_out =  np.asarray(symm_data[key])
        np.save(key+"_symm", symm_out)
        print key, symm_out.shape

def main(args):
    with open('Atom_Dict.json') as AtomsIn:
        AtomInputList = json.loads(AtomsIn.read())
    with open('Parameters.json') as params_in:
        params = json.loads(params_in.read())
    OutputDimension2 = generate_output_dimensions(params)
    keylist, All_G_Data, counter_dict = initialize_numpy_bins(AtomInputList, OutputDimension2)
    energy_dict = {}
    for wavefunction in keylist:
        print wavefunction
        labels, distance_matrix = retrieve_coordinates(wavefunction, params)
        #cycle through central atoms
        for i in range(0, len(distance_matrix)):
            if args.G1flag == True:
                G1_Data = []
            if args.G2flag == True:
                G2_Data = []
            if args.G3flag == True:
                G3_Data = []
            if args.G4flag == True:
                G4_Data = []
            if args.G5flag == True:
                G5_Data = []
            #retain atom label for central atom
            atom_type = filter(lambda x: not x.isdigit(), labels[i])
            #TODO: Pass around numpy array for operations instead of doing operations on single numbers
            if atom_type not in All_G_Data.keys():
                raise Exception(
                    'Unrecognized atom detected in database, please include it in the input dictionary using the -a or --atoms flag')
            #cycle through cutoff values
            for cutoff in params['cutoff']:
                fcRij_matrix = calculate_fcRij_matrix(distance_matrix, cutoff)
                if args.G1flag == True:
                    G1_Row = calculate_g1(fcRij_matrix[i])
                    G1_Data.append(np.sum(G1_Row))
                if args.G2flag == True:
                    for width in params['gausswidth']:
                        for distance in params['radial_distance']:
                            G2_Row = calculate_g2(
                                fcRij_matrix[i], distance_matrix[i], width, distance)
                            G2_Data.append(np.sum(G2_Row))
                if args.G3flag == True:
                    for period in params['period_length']:
                        G3_Row = calculate_g3(
                            fcRij_matrix[i], distance_matrix[i], period)
                        G3_Data.append(np.sum(G3_Row))
                if args.G4flag == True or args.G5flag == True:
                    for value in params['lambda_value']:
                        for resolution in params['angular_resolution']:
                            for width in params['gausswidth']:
                                G4_Sum = 0
                                G5_Sum = 0
                                for j in range(0, len(distance_matrix)):
                                    #This is here so we only take the upper triangle of the matrix
                                    for k in range(0, len(distance_matrix)):
                                        if i == j or j == k or i == k:
                                            pass
                                        elif distance_matrix[i,k] < 1e-7 or distance_matrix[i,k] >= cutoff or distance_matrix[i,j] < 1e-7 or distance_matrix[i,j] > cutoff:
                                            pass
                                        else:
                                            angle = calculate_angle(distance_matrix[i, j], distance_matrix[i, k], distance_matrix[j, k])
                                            if args.G4flag == True:
                                                G4_Sum +=  calculate_g4(
                                                    fcRij_matrix[i, j], fcRij_matrix[i, k], calculate_fcRij_matrix(distance_matrix[j, k], cutoff), distance_matrix[i, j], distance_matrix[i, k], distance_matrix[j, k], angle, value, resolution, width)
                                            if args.G5flag == True:
                                                G5_Sum += calculate_g5(
                                                    fcRij_matrix[i, j], fcRij_matrix[i, k], distance_matrix[i, j], distance_matrix[i, k], angle, value, resolution, width)
                                G4_Data.append(G4_Sum)
                                G5_Data.append(G5_Sum)
            G1_Data = np.asarray(G1_Data)
            G2_Data = np.asarray(G2_Data)
            G3_Data = np.asarray(G3_Data)
            G4_Data = np.asarray(G4_Data)
            G5_Data = np.asarray(G5_Data)
            G_Master_Data = np.concatenate((G1_Data, G2_Data, G3_Data, G4_Data, G5_Data), axis =0)
            All_G_Data[atom_type][counter_dict[atom_type]] = G_Master_Data
            counter_dict[atom_type] += 1
        for intfile in wavefunction:
            atom_label, energy = gen_energy_list(intfile)
            if atom_label not in energy_dict:
                energy_dict[atom_label] = []
            energy_dict[atom_label].append(energy)
    save_energy_data(energy_dict)
    save_g_values(All_G_Data)


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
