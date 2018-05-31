from __future__ import division
import os
import math
import random
import numpy as np

from scipy.spatial.distance import pdist, squareform

HEnergyList = []
CEnergyList = []
OEnergyList = []
HydrogenCounter = 0
CarbonCounter = 0
OxygenCounter = 0
CarbonWritten = 0
OxygenWritten = 0
HydrogenWritten = 0
OutputDimension2 = 0


def calculate_fcRij(distance):
    global Cutoff
    if distance <= 1e-7:
        pass
    if distance > cutoff:
        fcRij = 0
    else:
        fcRij = 0.5 * (math.cos(math.pi * distance / Cutoff) + 1)
    return fcRij

# Removed because it's redundant if I calculate fcRij
# def calculate_g1(fcRij):
#     G1 = 0
#     G1 += fcRij
#     eturn G1


def calculate_g2(fcRij, distance):
    G2 = 0
    global gausswidth
    G2 += math.exp(gausswidth * (distance) ** 2) * fcRij
    return G2


def calculate_g3(fcRij, distance):
    G3 = 0
    global period_length
    G3 += math.cos(distance * period_length) * fcRij
    return G3


def calculate_g4(fcRij, distance, angle):
#There are problems with the angle calculations
    G4 = 0
    global gausswidth, Lamba, angular_resolution
    G4 += (1 + Lamba * math.cos(angle)) ** angular_resolution * math.exp(-gausswidth * (distance12 ** 2 + distance13 ** 2 + distance23 ** 2)) * calculate_fcRij(distance12) * calculate_fcRij(distance13) * calculate_fcRij(distance23)
    G4 = G4 * 2 ** (1 - angular_resolution)
    return G4


def calculate_g5(fcRij, distance, angle):
#There are problems with the angle calculations
    G5 = 0
    global gausswidth, Lamba, angular_resolution
    G5 += (1 + Lamba * math.cos(angle)) ** angular_resolution * math.exp(-gausswidth * (distance12 ** 2 + distance13 ** 2)) * calculate_fcRij(distance12) * calculate_fcRij(distance13)
    G5 = G5 * 2 ** (1 - angular_resolution)
    return G5


def detect_angle(atomic_coordinates, atomic_coordinates_cutoff):
    neighbors_list = []
    neighbors_index_list = []
    for i, element in enumerate(atomic_coordinates_cutoff):
        if element >= 1e-7:
            neighbors_list.append(element)
            neighbors_index_list.append(i)
    if len(neighbors_list) == 2:
        distance12 = neighbors_list[0]
        distance13 = neighbors_list[1]
        distance23 = atomic_coordinates[neighbors_index_list[0], neighbors_index_list[1]]
        angle213 = calculate_angle(distance12, distance13, distance23)
        return (angle213,)
    if len(neighbors_list) == 3:
        distance12 = neighbors_list[0]
        distance13 = neighbors_list[1]
        distance14 = neighbors_list[2]
        distance23 = atomic_coordinates[neighbors_index_list[0], neighbors_index_list[1]]
        distance24 = atomic_coordinates[neighbors_index_list[0], neighbors_index_list[2]]
        distance34 = atomic_coordinates[neighbors_index_list[1], neighbors_index_list[2]]
        angle213 = calculate_angle(distance12, distance13, distance23)
        angle214 = calculate_angle(distance12, distance14, distance24)
        angle314 = calculate_angle(distance13, distance14, distance34)
        return (angle213, angle214, angle314)


def calculate_angle(distance12, distance13, distance23):
    angle = math.acos((distance12 ** 2 + distance13 ** 2 - distance23 ** 2) / (2 * distance12 * distance13))
    return angle


def gen_energy_list(int_file):
    with open(int_file, 'r') as atomic_file:
        atomic_lines = atomic_file.readlines()
        for line in atomic_lines:
            if line == atomic_lines[6]:
                atom_label = line[0]
            if line.startswith('              K'):
                floatenergy = float(line.split()[3])
    if atom_label == 'C':
        CEnergyList.append(floatenergy)
    if atom_label == 'H':
        HEnergyList.append(floatenergy)
    if atom_label == 'O':
        OEnergyList.append(floatenergy)


def label_to_charge(label):
    atom_label = label[0]
    if atom_label == 'H':
        return 1
    if atom_label == 'C':
        return 6
    if atom_label == 'O':
        return 8


def retrieve_coordinates(wavefunction, Cutoff):
    with open(wavefunction, 'r') as waveinput:
        wavelines = waveinput.readlines()
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
    matrix_cutoff[matrix_cutoff >= Cutoff] = 0
    return label_list, distance_matrix, matrix_cutoff

def gen_symm_functions(matrix, labels, matrix_cutoff):
    for i in range(0, len(matrix)):
        global CarbonWritten, OxygenWritten, HydrogenWritten
        data = []
        fcRij_list = []
        fcRij_total = 0.0
        for distance in matrix[i]:
            if distance > 1e-7:
                symm_function.append(distance)
                fcRij_list.append(calculate_fcRij(distance))
            atom_label = labels[i][0]
        for fcRij in fcRij_list:
            if G1flag == 'Y':
                G1_Total += fcRij
            if G2flag == 'Y':
                G2_Total += calculate_g2(fcRij, distance)
            if G3flag == 'Y':
                G3_Total += calculate_g3(fcRij, distance)
            if G4flag == 'Y':
                angles = detect_angle(matrix, matrix_cutoff)
                for angle in angles:
                    G4_Total += calculate_g4(fcRij, distance, angle)
            if G5flag == 'Y':
                angles = detect_angle(matrix, matrix_cutoff)
                for angle in angles:
                    G5_Total += calculate_g5(fcRij, distance, angle)
        if atom_label == 'C':
            CarbonArray[CarbonWritten] = symm_function
            CarbonWritten+=1
        if atom_label == 'H':
            HydrogenArray[HydrogenWritten] = symm_function
            HydrogenWritten+=1
        if atom_label == 'O':
            OxygenArray[OxygenWritten] = symm_function
            OxygenWritten += 1


Cutoff = input('Enter Cutoff in angstroms: ')
G1flag = input('Use G1? (Y/N)')
G2flag = input('Use G2? (Y/N)')
G3flag = input('Use G3? (Y/N)')
G4flag = input('Use G4? (Y/N)')
G5flag = input('Use G5? (Y/N)')
if G1flag:
    OutputDimension2 += 1
if G2flag:
    OutputDimension2 += 1
if G3flag:
    OutputDimension2 += 1
if G4flag:
    OutputDimension2 += 1
if G5flag:
    OutputDimension2 += 1

master_dict = {}

for atomfilename in filelist:
    if atomfilename.endswith('int'):
        wavefunction = atomfilename[0:21]
        print str(wavefunction)
        if wavefunction not in master_dict:
            atomfilelist = [atomfilename]
            master_dict[wavefunction] = atomfilelist
        else:
            atomfilelist = master_dict[wavefunction]
            atomfilelist.append(atomfilename)
            master_dict[wavefunction] = atomfilelist
        if atomfilename.startswith('dsC7O2H10nsd') and ".wfnC" in atomfilename and atomfilename.endswith('.int'):
            CarbonCounter += 1
        if atomfilename.startswith('dsC7O2H10nsd') and ".wfnH" in atomfilename and atomfilename.endswith('.int'):
            HydrogenCounter += 1
        if atomfilename.startswith('dsC7O2H10nsd') and ".wfnO" in atomfilename and atomfilename.endswith('.int'):
            OxygenCounter += 1

keylist = master_dict.keys()
keylist.sort()
for wavefunction in keylist:
    print wavefunction
    labels, distance_matrix = retrieve_coordinates(wavefunction)
    gen_symm_functions(distance_matrix, labels)
    for intfile in master_dict[wavefunction]:
        gen_energy_list(intfile)

HydrogenArray = np.zeros([HydrogenCounter, OutputDimension2])
CarbonArray = np.zeros([CarbonCounter, OutputDimension2])
OxygenArray = np.zeros([OxygenCounter, OutputDimension2])

HydrogenEnergyOut = np.asarray(HEnergyList)
CarbonEnergyOut = np.asarray(CEnergyList)
OxygenEnergyOut = np.asarray(OEnergyList)

print "Hydrogen", HydrogenArray.shape, HydrogenEnergyOut.shape
print "Carbon", CarbonArray.shape, CarbonEnergyOut.shape
print "Oxygen", OxygenArray.shape, OxygenEnergyOut.shape

np.save(('H_Out_Symmetry_Functions'), HydrogenArray)
np.save(('C_Out_Symmetry_Functions'), CarbonArray)
np.save(('O_Out_Symmetry_Functions'), OxygenArray)
np.save(('Energy_H_Out_Symmetry_Functions'), HydrogenEnergyOut)
np.save(('Energy_C_Out_Symmetry_Functions'), CarbonEnergyOut)
np.save(('Energy_O_Out_Symmetry_Functions'), OxygenEnergyOut)
