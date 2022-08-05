import math
import pandas as pd
import decimal
import numpy as np


def energies_scan(name, conv=True):
    key = "SCF Done:"
    key2 = "Optimized Parameters"
    key3 = "The following ModRedundant input section has been read:"
    file = open(name, "r")
    lines = file.readlines()
    app = []
    app2 = []
    csv = []
    file.close()
    for number, line in enumerate(lines, 0):
        if key in line:
            app.append(number)
    for number, line in enumerate(lines, 0):
        if key2 in line:
            app2.append(number)
    arr = np.array(app)
    ile2 = len(app2)

    # Appending list with energies for each of the steps
    for x in range(0, ile2):
        num = int(app2[x])
        num2 = arr[arr < num].max()
        eneline = lines[num2]
        ene = eneline[26:40]
        ene2 = ene + " \n"
        csv.append(ene2)

    if conv == True:
        header = 'Energy difference [kcal/mol]'
        csv2 = []
        for kk in range(0, len(csv)):
            csv2.append((float(csv[kk].replace("\n", ""))
                         - float(csv[0].replace("\n", ""))) * 627.51)

    elif conv == False:
        header = 'Energy [Hartree]'
        csv2 = []
        for kk in range(0, len(csv)):
            csv2.append(float(csv[kk].replace("\n", "")))

    for number, line in enumerate(lines, 0):
        if key3 in line:
            z = number + 2
            t = number + 1

    step = lines[z].replace(lines[z][0: lines[z].index('S') + 1], '')
    for y in step:
        if y == ' ':
            step = step.replace(y, '', 1)
        else:
            break
    step = step.replace(step[0: step.index(' ') + 1], '', 1)
    for y in step:
        if y  == ' ':
            step = step.replace(y , '', 1)
        else:
            break
    step = float(step)

    steps = []
    for x in range(0, len(csv2)):
        steps.append(round(step * x, 3))

    result = pd.DataFrame({
        'Bond Lenght': steps,
        header: csv2
    })

    # Reading bridge atoms
    list = []
    atom = lines[t]
    for y in atom:
        if y == ' ' or y == 'A':
            atom = atom.replace(y, '', 1)
        else:
            break
    list.append(int(atom[0: atom.index(' ')]))
    atom = atom.replace(str(list[0]), '', 1)
    for y in atom:
        if y == ' ':
            atom = atom.replace(y, '', 1)
        else:
            break
    list.append(int(atom[0: atom.index(' ')]))
    atom = atom.replace(str(list[1]), '',  1)
    for y in atom:
        if y == ' ':
            atom = atom.replace(y, '', 1)
        else:
            break
    list.append(int(atom[0: atom.index(' ')]))
    return result, list

def atomicnumber_to_name(atomic_number, nonspecified_num=None, nonspecified_lbl=None):
    """
        Function changing atomic number to atom type label

        Parameters
        ----------
        atomic_number : str
            str of characters representing atomic number
        nonspecified_num: int
            definable additional atomic number, for those not included in ths function
        nonspecified_lbl: int
            definable additional atom type for nonspecified_num atomic number
        Returns
        -------
        str
            Atom type
        """
    if atomic_number == "1":
        atomic_numberx = "H"
    elif atomic_number == "2":
        atomic_numberx = "He"
    elif atomic_number == "3":
        atomic_numberx = "Li"
    elif atomic_number == "4":
        atomic_numberx = "Be"
    elif atomic_number == "5":
        atomic_numberx = "B"
    elif atomic_number == "6":
        atomic_numberx = "C"
    elif atomic_number == "7":
        atomic_numberx = "N"
    elif atomic_number == "8":
        atomic_numberx = "O"
    elif atomic_number == "9":
        atomic_numberx = "F"
    elif atomic_number == "10":
        atomic_numberx = "Ne"
    elif atomic_number == "11":
        atomic_numberx = "Na"
    elif atomic_number == "12":
        atomic_numberx = "Mg"
    elif atomic_number == "13":
        atomic_numberx = "Al"
    elif atomic_number == "14":
        atomic_numberx = "Si"
    elif atomic_number == "15":
        atomic_numberx = "P"
    elif atomic_number == "16":
        atomic_numberx = "S"
    elif atomic_number == "17":
        atomic_numberx = "Cl"
    elif atomic_number == "18":
        atomic_numberx = "Ar"
    elif atomic_number == "19":
        atomic_numberx = "K"
    elif atomic_number == "20":
        atomic_numberx = "Ca"
    elif atomic_number == "26":
        atomic_numberx = "Fe"
    elif atomic_number == "28":
        atomic_numberx = "Ni"
    elif atomic_number == "29":
        atomic_numberx = "Cu"
    elif atomic_number == "25":
        atomic_numberx = "Mn"
    elif atomic_number == "78":
        atomic_numberx = "Pt"
    elif atomic_number == "35":
        atomic_numberx = "Br"
    elif atomic_number == "53":
        atomic_numberx = "I"
    elif atomic_number == "47":
        atomic_numberx = "Ag"
    elif atomic_number == "47":
        atomic_numberx = "Au"
    elif atomic_number == str(nonspecified_num):
        atomic_numberx = str(nonspecified_lbl)

    return atomic_numberx


def extraction_atom_dataframe(name, file_type):
    """
        Function extracting atom matrix with 4 columns
        (Atom type, X, Y, Z coordinates)

        Parameters
        ----------
        name : str
            filename used to extract from
        file_type: str
            'opt' for logiles from optizmiation runs
            'scan' for logiles from scan runs
        Returns
        -------
        distance_frame: pd.DataFrame
            DataFrame containing all distances calculated by Gaussian, in row format:
            Atom 1, Atom 2, Distance
        angle_frame: pd.DataFrame
            DataFrame containing all angles calculated by Gaussian, in row format:
            Atom 1, Atom 2, Atom 3, Distance
        """
    if file_type == "scan":
        file = open(name, "r")
        lines = file.readlines()
        file.close()
        key = 'Symbolic Z-matrix:'
        key2 = "The following ModRedundant input section has been read:"
        for number, line in enumerate(lines, 0):  # Line before atom-type matrix
            if key in line:
                a_matrix_start = int(number) + 1
        for number, line in enumerate(lines, 0):  # Line after atom-type matrix
            if key2 in line:
                a_matrix_end = int(number) - 1
                break
        atom_matrix = lines[a_matrix_start:a_matrix_end]  # Creating atoms matrix
        atom_matrix[0] = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"

        '''Creating DataFrame of atoms and positions using atom_matrix.
        Each line is split into distinct elements, which are appended into
        related columns. Then, Dataframe is built.'''
        index = []
        col1 = []
        col2 = []
        col3 = []
        col4 = []
        for x in range(1, len(atom_matrix)):
            index.append(x)
            atom_matrix[x] = atom_matrix[x].replace(' ', '', 1)
            col1.append(atom_matrix[x][0: atom_matrix[x].index(' ')])
            atom_matrix[x] = atom_matrix[x].replace(col1[x - 1], '', 1)
            for y in atom_matrix[x]:
                if y == ' ':
                    atom_matrix[x] = atom_matrix[x].replace(' ', '', 1)
                else:
                    break
            col2.append(atom_matrix[x][0: atom_matrix[x].index(' ')])
            atom_matrix[x] = atom_matrix[x].replace(col2[x - 1], '', 1)
            for y in atom_matrix[x]:
                if y == ' ':
                    atom_matrix[x] = atom_matrix[x].replace(' ', '', 1)
                else:
                    break
            col3.append(atom_matrix[x][0: atom_matrix[x].index(' ')])
            atom_matrix[x] = atom_matrix[x].replace(col3[x - 1], '', 1)
            for y in atom_matrix[x]:
                if y == ' ':
                    atom_matrix[x] = atom_matrix[x].replace(' ', '', 1)
                else:
                    break
            col4.append(atom_matrix[x][0: atom_matrix[x].index(' ')])

        for x in range(0, len(col2)):
            col2[x] = float(col2[x])
            col3[x] = float(col3[x])
            col4[x] = float(col4[x])

        atom_dataframe = pd.DataFrame({
            'Atom type': col1,
            'x': col2,
            'y': col3,
            'z': col4
        }, index=index)
        return atom_dataframe

    elif file_type == 'opt':
        key = "Center     Atomic      Atomic             Coordinates (Angstroms)"
        key2 = "symmetry adapted cartesian basis functions of"
        file = open(name, "r")
        lines = file.readlines()
        file.close()
        app = []
        app2 = []

        for number, line in enumerate(lines, 0):
            if key in line:
                app.append(number)
        for number, line in enumerate(lines, 0):
            if key2 in line:
                app2.append(number)
        # Choosing lines for the last, optimized structure
        final1 = app[-1]
        final2 = app2[-1]

        # Building atom matrix
        start = int(final1) + 3
        finish = int(final2) - 3
        atom_matrix = lines[start:finish]

        '''Creating DataFrame of atoms and positions using atom_matrix.
        Each line is split into distinct elements, which are appended into
        related columns. Then, Dataframe is built.'''
        rowindex = []
        col1 = []
        col2 = []
        col3 = []
        col4 = []
        col_trash = []
        for x in range(0, len(atom_matrix)):
            for y in atom_matrix[x]:
                if y == ' ':
                    atom_matrix[x] = atom_matrix[x].replace(' ', '', 1)
                else:
                    break

            rowindex.append(int(atom_matrix[x][0: atom_matrix[x].index(' ')]))
            atom_matrix[x] = atom_matrix[x].replace(str(rowindex[x]), '', 1)

            for y in atom_matrix[x]:
                if y == ' ':
                    atom_matrix[x] = atom_matrix[x].replace(' ', '', 1)
                else:
                    break

            col1.append(int(atom_matrix[x][0: atom_matrix[x].index(' ')]))
            atom_matrix[x] = atom_matrix[x].replace(str(col1[x]), '', 1)

            for y in atom_matrix[x]:
                if y == ' ':
                    atom_matrix[x] = atom_matrix[x].replace(' ', '', 1)
                else:
                    break

            col_trash.append(int(atom_matrix[x][0: atom_matrix[x].index(' ')]))
            atom_matrix[x] = atom_matrix[x].replace(str(col_trash[x]), '', 1)

            for y in atom_matrix[x]:
                if y == ' ':
                    atom_matrix[x] = atom_matrix[x].replace(' ', '', 1)
                else:
                    break
            linecheck = atom_matrix[x]
            col2.append(atom_matrix[x][0: atom_matrix[x].index(' ')])
            atom_matrix[x] = atom_matrix[x].replace(str(col2[x]), '', 1)

            for y in atom_matrix[x]:
                if y == ' ':
                    atom_matrix[x] = atom_matrix[x].replace(' ', '', 1)
                else:
                    break

            col3.append(atom_matrix[x][0: atom_matrix[x].index(' ')])
            atom_matrix[x] = atom_matrix[x].replace(str(col3[x]), '', 1)

            for y in atom_matrix[x]:
                if y == ' ':
                    atom_matrix[x] = atom_matrix[x].replace(' ', '', 1)
                else:
                    break

            col4.append(float(atom_matrix[x]))

        for x in range(0, len(col1)):
            col1[x] = atomicnumber_to_name(str(col1[x]))
        for x in range(0, len(col2)):
            col2[x] = float(col2[x])
            col3[x] = float(col3[x])
            col4[x] = float(col4[x])

        atom_dataframe = pd.DataFrame({
            'Atom type': col1,
            'x': col2,
            'y': col3,
            'z': col4
        }, index=rowindex)
        return atom_dataframe

    else:
        print('Wrong type!')


def extraction_parameters(name, filetype):
    """Function extracting distance matrix with 3 columns
    (Atom 1, Atom 2, Value)
    as well as angle matrix with 4 columns
    (Atom 1, Atom 2, Atom 2, Value)

    Parameters
    ----------
    name : str
        Filename used to extract from
    Returns
    -------
    distance_frame: pd.DataFrame
        DataFrame containing all distances calculated by Gaussian, in row format:
        Atom 1, Atom 2, Distance
    angle_frame: pd.DataFrame
        DataFrame containing all angles calculated by Gaussian, in row format:
        Atom 1, Atom 2, Atom 3, Distance
    """

    file = open(name, "r")
    lines = file.readlines()
    file.close()
    if filetype == 'scan':
        key = "Initial Parameters"
        key2 = 'Trust Radius='
        for number, line in enumerate(lines, 0):  # Line after atom-type matrix
            if key in line:
                p_matrix_start = int(number) + 5
                break
        for number, line in enumerate(lines, 0):  # Finding line after parameter matrix
            if key2 in line:
                p_matrix_end = int(number) - 1
                break
        parameter_matrix = lines[p_matrix_start:p_matrix_end]
    if filetype == 'opt':
        app = []
        key = "Optimized Parameters"
        key2 = 'Non-Optimized Parameters'
        for number, line in enumerate(lines, 0):  # Line after atom-type matrix
            if key in line and key2 not in line:
                app.append(int(number) + 5)
                break
        p_matrix_start = app[-1]
        for x in range(app[-1], len(lines)):
            if 'A' not in lines[x][0: 10] and 'D' not in lines[x][0: 10] and 'R' not in lines[x][0: 10]:
                p_matrix_end = x
                break
        parameter_matrix = lines[p_matrix_start:p_matrix_end]
    rowindex_dist = []
    rowindex_ang = []
    rowindex_dih = []
    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    col6 = []
    col7 = []
    col8 = []
    col9 = []
    col10 = []
    col11 = []
    col12 = []
    for x in range(0, len(parameter_matrix)):
        if 'R' in parameter_matrix[x][0:10]:
            for y in parameter_matrix[x]:
                if y == ' ' or y == '!' or y == 'R':
                    parameter_matrix[x] = parameter_matrix[x].replace(y, '', 1)
                else:
                    break
            rowindex_dist.append(int(parameter_matrix[x][0: parameter_matrix[x].index(' ')]))
            parameter_matrix[x] = parameter_matrix[x].replace(str(rowindex_dist[x]), '', 1)

            for y in parameter_matrix[x]:
                if y == ' ' or y == 'R' or y == '(':
                    parameter_matrix[x] = parameter_matrix[x].replace(y, '', 1)
                else:
                    break
            col1.append(int(parameter_matrix[x][0: parameter_matrix[x].index(',')]))
            parameter_matrix[x] = parameter_matrix[x].replace(str(col1[x]) + ',', '', 1)

            col2.append(int(parameter_matrix[x][0: parameter_matrix[x].index(')')]))
            parameter_matrix[x] = parameter_matrix[x].replace(str(col2[x]) + ')', '', 1)
            for y in parameter_matrix[x]:
                if y == ' ':
                    parameter_matrix[x] = parameter_matrix[x].replace(' ', '', 1)
                else:
                    break
            col3.append(float(parameter_matrix[x][0: parameter_matrix[x].index(' ')]))

        elif 'A' in parameter_matrix[x][0:10]:
            for y in parameter_matrix[x]:
                if y == ' ' or y == '!' or y == 'A':
                    parameter_matrix[x] = parameter_matrix[x].replace(y, '', 1)
                else:
                    break
            rowindex_ang.append(int(parameter_matrix[x][0: parameter_matrix[x].index(' ')]))
            parameter_matrix[x] = parameter_matrix[x].replace(str(rowindex_ang[x - len(col1)]), '', 1)

            for y in parameter_matrix[x]:
                if y == ' ' or y == 'A' or y == '(':
                    parameter_matrix[x] = parameter_matrix[x].replace(y, '', 1)
                else:
                    break
            col4.append(int(parameter_matrix[x][0: parameter_matrix[x].index(',')]))
            parameter_matrix[x] = parameter_matrix[x].replace(str(col4[x - len(col1)]) + ',', '', 1)
            col5.append(int(parameter_matrix[x][0: parameter_matrix[x].index(',')]))
            parameter_matrix[x] = parameter_matrix[x].replace(str(col5[x - len(col1)]) + ',', '', 1)


            col6.append(int(parameter_matrix[x][0: parameter_matrix[x].index(')')]))
            parameter_matrix[x] = parameter_matrix[x].replace(str(col6[x - len(col1)]) + ')', '', 1)
            for y in parameter_matrix[x]:
                if y == ' ':
                    parameter_matrix[x] = parameter_matrix[x].replace(' ', '', 1)
                else:
                    break
            col7.append(float(parameter_matrix[x][0: parameter_matrix[x].index(' ')]))

        elif 'D' in parameter_matrix[x][0:10]:
            for y in parameter_matrix[x]:
                if y == ' ' or y == '!' or y == 'D':
                    parameter_matrix[x] = parameter_matrix[x].replace(y, '', 1)
                else:
                    break
            line_help = parameter_matrix[x]
            rowindex_dih.append(int(parameter_matrix[x][0: parameter_matrix[x].index(' ')]))
            parameter_matrix[x] = parameter_matrix[x].replace(str(rowindex_dih[x - (len(col1) + len(col4))]), '', 1)

            for y in parameter_matrix[x]:
                if y == ' ' or y == 'D' or y == '(':
                    parameter_matrix[x] = parameter_matrix[x].replace(y, '', 1)
                else:
                    break
            col8.append(int(parameter_matrix[x][0: parameter_matrix[x].index(',')]))
            parameter_matrix[x] = parameter_matrix[x].replace(str(col8[x - len(col1) - len(col4)]) + ',', '', 1)
            col9.append(int(parameter_matrix[x][0: parameter_matrix[x].index(',')]))
            parameter_matrix[x] = parameter_matrix[x].replace(str(col9[x - len(col1) - len(col4)]) + ',', '', 1)
            col10.append(int(parameter_matrix[x][0: parameter_matrix[x].index(',')]))
            parameter_matrix[x] = parameter_matrix[x].replace(str(col10[x - len(col1) - len(col4)]) + ',', '', 1)

            col11.append(int(parameter_matrix[x][0: parameter_matrix[x].index(')')]))
            parameter_matrix[x] = parameter_matrix[x].replace(str(col11[x - len(col1) - len(col4)]) + ')', '', 1)
            for y in parameter_matrix[x]:
                if y == ' ':
                    parameter_matrix[x] = parameter_matrix[x].replace(' ', '', 1)
                else:
                    break
            col12.append(float(parameter_matrix[x][0: parameter_matrix[x].index(' ')]))

    distance_frame = pd.DataFrame({
        'Atom 1': col1,
        'Atom 2': col2,
        'Distance': col3
    }, index=rowindex_dist)
    angle_frame = pd.DataFrame({
        'Atom 1': col4,
        'Atom 2': col5,
        'Atom 3': col6,
        'Angle': col7,
    }, index=rowindex_ang)
    dihedral_frame = pd.DataFrame({
        'Atom 1': col8,
        'Atom 2': col9,
        'Atom 3': col10,
        'Atom 4': col11,
        'Dihedral': col12,
    }, index=rowindex_dih)
    return distance_frame, angle_frame, dihedral_frame

def float2string(f):
    ctx = decimal.Context()
    ctx.prec = 7
    d1 = ctx.create_decimal(repr(f))
    return format(d1, 'f')


def distance_calc(where, atom1_number, atom2_number):
    """
        Function returning distance between two points in XYZ defined space.

        Parameters
        ----------
        where: pd.DataFrame
            DataFrame with XYZ coordinates in columns, given as 'x', 'y', 'z'
        atom1_number : int
           Number of point 1
        atom2_number: int
           Number of point 2

        Returns
        -------
        distance: float
       """

    distance = math.sqrt(((where.loc[atom1_number, 'x'] - where.loc[atom2_number, 'x']) ** 2)
              + ((where.loc[atom1_number, 'y'] - where.loc[atom2_number, 'y']) ** 2)
              + ((where.loc[atom1_number, 'z'] - where.loc[atom2_number, 'z']) ** 2))

    return distance

def full_matrix(name, filetype='opt'):
    """
        Function extracting 5 DataFrames from file, that can be further used
        in the other functions from GaussianMaster.

        Parameters
        ----------
        name : str
            name or path of a file to extract from
        filetype: str
            Type of logfile - if result of scan run then 'scan',
            if optimalization - 'opt' (default)

        Returns
        -------
        atom_matrix: pd.DataFrame
            Table with Atom, X, Y, Z in columns
        distance_frame: pd.DataFrame
            Table with Atom1, Atom2, Distance in columns (Atoms as numbers)
        distance_frame_clean: pd.DataFrame
            Table with Atom1, Atom2, Distance in columns (Atoms as types)
        angle_frame: pd.DataFrame
            Table with Atom1, Atom2, Atom3, Angle in columns (Atoms as numbers)
        angle_frame_clean: pd.DataFrame
            Table with Atom1, Atom2, Atom3, Angle in columns (Atoms as types)
        """
    if filetype != 'opt' and filetype != 'scan':
        raise NameError('Wrong filetype parameter')
    elif '.log' not in name:
        raise NameError('Wrong type of file! Only use .log files.')
    else:

        # Using submodules to extract number-based tables
        atom_matrix = extraction_atom_dataframe(name, filetype)
        distance_frame, angle_frame, dihedral_frame = extraction_parameters(name, filetype)

        # Cleaning data by turning numbers to their respective atom types
        distance_frame_clean = distance_frame.copy()
        angle_frame_clean = angle_frame.copy()
        dihedral_frame_clean = dihedral_frame.copy()
        r, c = distance_frame.shape
        for x in range(0, r):
            for y in range(0, 2):
                distance_frame_clean.iloc[x, y] = atom_matrix.iloc[int(distance_frame.iloc[x, y]) - 1, 0]
        r, c = angle_frame.shape
        for x in range(0, r):
            for y in range(0, 3):
                angle_frame_clean.iloc[x, y] = atom_matrix.iloc[int(angle_frame.iloc[x, y]) - 1, 0]
        r, c = dihedral_frame.shape
        for x in range(0, r):
            for y in range(0, 4):
                dihedral_frame_clean.iloc[x, y] = atom_matrix.iloc[int(dihedral_frame.iloc[x, y]) - 1, 0]

        return atom_matrix, distance_frame, distance_frame_clean, angle_frame, angle_frame_clean, dihedral_frame, dihedral_frame_clean


def bridge_atom_selection(atom_matrix, distance_matrix, distance_matrix_clean, donor='O', acceptor='N'):
    """
        Function selecting three atoms composing a hydrogen bridges.
        It returns three lists containing atom type (donor, proton, acceptor), 
        with corresponding elements of the lists being on the same place in each one.
        

        Parameters
        ----------
        atom_matrix : pd.DataFrame
             atom and its xyz coordinates matrix from full_matrix()
        distance_matrix: pd.DataFrame
             atom1 number, atom2 number and distance value matrix from full_matrix()
        distance_matrix_clean: pd.DataFrame
             atom1 type, atom2 type and distance value matrix from full_matrix()
        donor: string
            type of donor atom, default O
        acceptor: string
            type of donor atom, default N

        Returns
        -------
        donor_atoms: list
            List of donor atom numbers.
        protons: list
            List of proton numbers.
        acceptor_atoms: list
            List of acceptor atom numbers.
        """

    donor_atoms = []
    protons = []
    acceptor_atoms = []
    for x in range(0, len(distance_matrix_clean)):
        if distance_matrix_clean.iloc[x, 0] == donor and distance_matrix_clean.iloc[x, 1] == 'H':
            donor_atoms.append(distance_matrix.iloc[x, 0])
            protons.append(distance_matrix.iloc[x, 1])
        elif distance_matrix_clean.iloc[x, 0] == 'H' and distance_matrix_clean.iloc[x, 1] == donor:
            donor_atoms.append(distance_matrix.iloc[x, 1])
            protons.append(distance_matrix.iloc[x, 0])

    # Checking if acceptor of a given type in range
    for x in range(0, len(protons)):
        for y in range(1, atom_matrix.index[-1]):
            if y == protons[x] or y in donor_atoms:
                pass
            elif atom_matrix.loc[y, 'Atom type'] == acceptor:
                if distance_calc(atom_matrix, protons[x], y)  <= 2.5:
                    acceptor_atoms.append(atom_matrix.index[y - 1])

    return donor_atoms, protons, acceptor_atoms


def bridge_parameters(donor_atoms, protons, acceptor_atoms, atom_matrix, which=1):
    """
        Function selecting three atoms composing a single, specified hydrogen bridge.

        Parameters
        ----------
        donor_atoms : list
             List containing donor atoms
        protons: list
             List containing protons
        acceptor_atoms: list
             List containing acceptor atoms
        atom_matrix: pd.DataFrame
            Atom matrix from full_matrix()
        which: int
            Used in case of multiple potential birdge protons, 1 by default.

        Returns
        -------
        acceptor_atom: int
            Acceptor atom number for a given bridge.
        proton: int
            Proton atom number for a given bridge.
        donor_atom: int
            Donor atom number for a given bridge.
        """

    acceptor_atom = 'None'
    donor_atom = 'None'
    proton = protons[which - 1]
    # Checking distance for proton to acceptor and donor
    for x in range(0, len(acceptor_atoms)):
        if distance_calc(atom_matrix, protons[which - 1], acceptor_atoms[x]) <= 2.5:
            acceptor_atom = acceptor_atoms[x]
            break
    for x in range(len(donor_atoms)):
        if distance_calc(atom_matrix, protons[which - 1], donor_atoms[x]) <= 2.5 \
                and str(donor_atoms[x]) != str(acceptor_atom):
            donor_atom = donor_atoms[x]

    return acceptor_atom, proton, donor_atom

def com_builder(what, functional_plus_basis, comtype, name, steps=20, dist=0.05, a_bridge=0, p_bridge=0, d_bridge=0, charge=0, multiplicity=1):
    """
            Function preparing xyz string to be written as .com file

            Parameters
            ----------
            what : pd.DataFrame
                List containing donor atoms
            functional_plus_basis: str
                String in the format '# func/basis,
                i.e. # PBEPBE/6-311++G(2d,2p)
            comtype: str
                String with com type that will be the result, changing xyz contents:
                'scan' - com for scans
                'opt' - com for optimization with frequencies
                'wfn' - combined wfn + Hirshfeld populations
            name: str
                Name to be written into comfile
            steps: int
                (optional) number of steps for scan, by default 20
            dist: int
                (optional) step distance for scan, by default 0.05
            a_bridge: int
                (optional) number of acceptor atom for scans
            p_bridge: int
                (optional) number of proton atom for scans
            d_bridge: int
                (optional) number of donor atom for scans
            charge: int
                (optional) charge, by default 0
            charge: int
                (optional) multiplicity, by default 0

            Returns
            -------
            xyz: str
                String ready to be written as .com
            """
    if comtype == 'scan':
        xyz = functional_plus_basis + ' opt=modredundant \n \n' \
              + str(name) + '\n \n' + str(charge) + ' ' + str(multiplicity) + '\n'
        for x in range(0, len(what)):
            xyz += str(what.iloc[x, 0])
            help_string = '               ' + float2string(what.iloc[x, 1])
            help_string = help_string.replace(' ', '', len(str(what.iloc[x, 0]))
                                              + float2string(what.iloc[x, 1]).count('-'))
            xyz += help_string
            help_string = '           ' + float2string(what.iloc[x, 2])
            help_string = help_string.replace(' ', '',
                                              (len(float2string(what.iloc[x, 1]))
                                               - float2string(what.iloc[x, 1]).count('-'))
                                              + float2string(what.iloc[x, 2]).count('-'))
            xyz += help_string
            help_string = '           ' + float2string(what.iloc[x, 3])
            help_string = help_string.replace(' ', '',
                                              (len(float2string(what.iloc[x, 2]))
                                               - float2string(what.iloc[x, 2]).count('-'))
                                              + float2string(what.iloc[x, 3]).count('-'))
            xyz += help_string + '\n'
        xyz += '\n' +  str(d_bridge) + ' ' +  str(p_bridge) + ' ' +  str(a_bridge) + ' F\n'
        xyz +=  str(d_bridge) + ' ' +  str(p_bridge) + ' S ' + str(steps) + ' ' + str(dist) + '\n\n\n'
        return xyz

    elif comtype == 'opt':
        xyz = functional_plus_basis + ' opt freq \n \n' \
              + str(name) + '\n \n' + str(charge) + ' ' + str(multiplicity) + '\n'
        for x in range(0, len(what)):
            xyz += str(what.iloc[x, 0])
            help_string = '               ' + float2string(what.iloc[x, 1])
            help_string = help_string.replace(' ', '',
                                              len(str(what.iloc[x, 0]))
                                              + float2string(what.iloc[x, 1]).count('-'))
            xyz += help_string
            help_string = '           ' + float2string(what.iloc[x, 2])
            help_string = help_string.replace(' ', '',
                                              (len(float2string(what.iloc[x, 1]))
                                               - float2string(what.iloc[x, 1]).count('-'))
                                              + float2string(what.iloc[x, 2]).count('-'))
            xyz += help_string
            help_string = '           ' + float2string(what.iloc[x, 3])
            help_string = help_string.replace(' ', '',
                                              (len(float2string(what.iloc[x, 2]))
                                               - float2string(what.iloc[x, 2]).count('-'))
                                              + float2string(what.iloc[x, 3]).count('-'))
            xyz += help_string + '\n'
        xyz += '\n\n\n'
        return xyz

    elif comtype == 'wfn':
        xyz = functional_plus_basis + ' sp scf=tight pop=hirshfeld out=wfn \n \n' \
              + str(name) + '\n \n' + str(charge) + ' ' + str(multiplicity) + '\n'
        for x in range(0, len(what)):
            xyz += str(what.iloc[x, 0])
            help_string = '               ' + float2string(what.iloc[x, 1])
            help_string = help_string.replace(' ', '', len(str(what.iloc[x, 0]))
                                              + float2string(what.iloc[x, 1]).count('-'))
            xyz += help_string
            help_string = '           ' + float2string(what.iloc[x, 2])
            help_string = help_string.replace(' ', '', (len(float2string(what.iloc[x, 1]))
                                                        - float2string(what.iloc[x, 1]).count('-'))
                                              + float2string(what.iloc[x, 2]).count('-'))
            xyz += help_string
            help_string = '           ' + float2string(what.iloc[x, 3])
            help_string = help_string.replace(' ', '', (len(float2string(what.iloc[x, 2]))
                                                        - float2string(what.iloc[x, 2]).count('-'))
                                              + float2string(what.iloc[x, 3]).count('-'))
            xyz += help_string + '\n'
        xyz += '\n' + name + '.wfn\n\n'
        return xyz
    else:
        raise NameError('Wrong comptype parameter')
