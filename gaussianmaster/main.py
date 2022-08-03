import math

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

    import gaussianmaster.core.extractor as ma

    # Using submodules to extract number-based tables
    atom_matrix = ma.extraction_atom_dataframe(name, filetype)
    distance_frame, angle_frame = ma.extraction_parameters(name)

    # Cleaning data by turning numbers to their respective atom types
    distance_frame_clean = distance_frame.copy()
    angle_frame_clean = angle_frame.copy()
    r, c = distance_frame.shape
    for x in range(0, r):
        for y in range(0, 2):
            distance_frame_clean.iloc[x, y] = atom_matrix.iloc[int(distance_frame.iloc[x, y]) - 1, 0]
    r, c = angle_frame.shape
    for x in range(0, r):
        for y in range(0, 3):
            angle_frame_clean.iloc[x, y] = atom_matrix.iloc[int(angle_frame.iloc[x, y]) - 1, 0]

    return atom_matrix, distance_frame, distance_frame_clean, angle_frame, angle_frame_clean


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
                if math.sqrt((atom_matrix.loc[protons[x], 'x'] - atom_matrix.loc[y, 'x']) ** 2
                             + (atom_matrix.loc[protons[x], 'y'] - atom_matrix.loc[y, 'y']) ** 2
                             + (atom_matrix.loc[protons[x], 'z'] - atom_matrix.loc[y, 'z']) ** 2) <= 2.5:
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
    # TODO Naprawic zjebane liczenie wektorow
    # Checking distance for proton to acceptor and donor
    for x in range(0, len(acceptor_atoms)):
        if math.sqrt(((atom_matrix.loc[protons[which - 1], 'x'] - atom_matrix.loc[acceptor_atoms[x] , 'x']) ** 2)
                     + ((atom_matrix.loc[protons[which - 1], 'y'] - atom_matrix.loc[acceptor_atoms[x], 'y']) ** 2)
                     + ((atom_matrix.loc[protons[which - 1], 'z'] - atom_matrix.loc[acceptor_atoms[x], 'z']) ** 2)) <= 2.5:
            acceptor_atom = acceptor_atoms[x]
            break
    for x in range(len(donor_atoms)):
        if math.sqrt(((atom_matrix.loc[protons[which - 1], 'x'] - atom_matrix.loc[donor_atoms[x] , 'x']) ** 2)
                     + ((atom_matrix.loc[protons[which - 1], 'y'] - atom_matrix.loc[donor_atoms[x], 'y']) ** 2)
                     + ((atom_matrix.loc[protons[which - 1], 'z'] - atom_matrix.loc[donor_atoms[x], 'z']) ** 2)
                     ) <= 2.5 and str(
            donor_atoms[x]) != str(acceptor_atom):

            donor_atom = donor_atoms[x]

    return acceptor_atom, proton, donor_atom




print('a')
