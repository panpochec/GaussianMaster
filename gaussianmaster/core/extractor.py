import pandas as pd


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


def extraction_parameters(name):
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
    key = "Initial Parameters"
    key2 = 'Trust Radius='
    for number, line in enumerate(lines, 0):  # Line after atom-type matrix
        if key in line:
            p_matrix_start = int(number) + 5
            break
    for number, line in enumerate(lines, 0):  # Finding line after parameter matrix
        if key2 in line:
            p_matrix_end = int(number) - 2
            break
    parameter_matrix = lines[p_matrix_start:p_matrix_end]

    rowindex_dist = []
    rowindex_ang = []
    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    col6 = []
    col7 = []
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

    return distance_frame, angle_frame

#atom_matrix = extraction_atom_dataframe('meta.log', 'opt')
