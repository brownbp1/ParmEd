"""
This package contains classes responsible for loading and dumping PyRosetta
Pose objects.
"""

from __future__ import print_function

from parmed.exceptions import RosettaError
from parmed.periodic_table import AtomicNum, Mass
from parmed.structure import Structure
from parmed.topologyobjects import Atom, ExtraPoint, Bond
from parmed.utils.six.moves import range

try:
    from pyrosetta import Pose, AtomID, rosetta
except ImportError:
    Pose = AtomID = rosetta = None


# Typedef for convenience
Vec3D = rosetta.numeric.xyzVector_double_t


# Convert residue atom id to raw atom count id
def _n_prior(pose, nbr):
    prior = -1
    for i in range(1, nbr.rsd()):
        prior += pose.residue(i).natoms()
    return prior + nbr.atomno()


#
def _map_cap_atoms_struct_to_pose(atom_name, res_name):
    if atom_name == "CH3" and res_name == "NME":
        return "CN"
    elif atom_name == "N" and res_name == "NME":
        return "NM"
    elif atom_name == "H" and res_name == "NME":
        return "HM"
    elif atom_name == "HH31" and res_name == "NME":
        return "1HN"
    elif atom_name == "HH32" and res_name == "NME":
        return "2HN"
    elif atom_name == "HH33" and res_name == "NME":
        return "3HN"
    elif atom_name == "C" and res_name == "ACE":
        return "CO"
    elif atom_name == "O" and res_name == "ACE":
        return "OP1"
    elif atom_name == "CH3" and res_name == "ACE":
        return "CP2"
    elif atom_name == "H1" and res_name == "ACE":
        return "1HP2"
    elif atom_name == "H2" and res_name == "ACE":
        return "2HP2"
    elif atom_name == "H3" and res_name == "ACE":
        return "3HP2"
    return None


def _map_cap_atoms_pose_to_struct(atom_name, res_name):
    if atom_name == "CN" and res_name == "ACE":
        return "CH3"
    elif atom_name == "NM" and res_name == "ACE":
        return "N"
    elif atom_name == "CO" and res_name == "NME":
        return "C"
    elif atom_name == "OP1" and res_name == "NME":
        return "O"
    elif atom_name == "CP2" and res_name == "NME":
        return "CH3"
    return None


#
def _map_atoms_struct_to_pose(key, res_name=None):

    # TODO: this should really be a generalized BCL substructure match at the residue level in C++ Rosetta
    # to create maps on-the-fly;
    # Amber --> Rosetta backbone atom map

    # Special cases for THR and aromatic ring residues
    if res_name == "THR" and key == "HG1":
        return "HG1"
    elif (res_name == "PHE" or res_name == "TYR" or res_name == "HIE" or res_name == "HID") and key == "HD1":
        return "HD1"
    elif (res_name == "PHE" or res_name == "TYR" or res_name == "HIE" or res_name == "HID" or res_name == "ASH") and key == "HD2":
        return "HD2"
    elif (res_name == "PHE" or res_name == "TYR" or res_name == "HIE" or res_name == "HID") and key == "HE1":
        return "HE1"
    elif (res_name == "PHE" or res_name == "TYR" or res_name == "HIE" or res_name == "HID") and key == "HE2":
        return "HE2"
    elif res_name == "TRP" and key == "HD1":
        return "HD1"
    elif res_name == "TRP" and key == "HE1":
        return "HE1"
    elif res_name == "TRP" and key == "HE3":
        return "HE3"
    elif res_name == "TRP" and key == "HZ2":
        return "HZ2"
    elif res_name == "TRP" and key == "HZ3":
        return "HZ3"
    elif res_name == "TRP" and key == "HH2":
        return "HH2"

    # Generic
    amber_to_rosetta_atom_names = {
        'N': 'N',
        'H': 'H',
        'CA': 'CA',
        'C': 'C',
        'O': 'O',

        'HA2': '1HA',
        'HA3': '2HA',

        'HB1': '3HB',
        'HB2': '1HB',
        'HB3': '2HB',

        'HG1': '3HG',
        'HG1': '3HG',
        'HG2': '1HG',
        'HG3': '2HG',

        'HD1': '3HD',
        'HD2': '1HD',
        'HD3': '2HD',

        'HE1': '3HE',
        'HE2': '1HE',
        'HE3': '2HE',

        'HH1': '3HH',
        'HH2': '1HH',
        'HH3': '2HH',

        'HZ1': '3HZ',
        'HZ2': '1HZ',
        'HZ3': '2HZ',

        # Branched sidechains
        'HB11': '1HB1',
        'HB12': '2HB1',
        'HB13': '3HB1',
        'HB21': '1HB2',
        'HB22': '2HB2',
        'HB23': '3HB2',

        'HG11': '3HG1',
        'HG12': '1HG1',
        'HG13': '2HG1',
        'HG21': '1HG2',
        'HG22': '2HG2',
        'HG23': '3HG2',

        'HD11': '1HD1',
        'HD12': '2HD1',
        'HD13': '3HD1',
        'HD21': '1HD2',
        'HD22': '2HD2',
        'HD23': '3HD2',

        'HE11': '1HE1',
        'HE12': '2HE1',
        'HE13': '3HE1',
        'HE21': '1HE2',
        'HE22': '2HE2',
        'HE23': '3HE2',

        'HH11': '1HH1',
        'HH12': '2HH1',
        'HH13': '3HH1',
        'HH21': '1HH2',
        'HH22': '2HH2',
        'HH23': '3HH2',

        'HZ11': '1HZ1',
        'HZ12': '2HZ1',
        'HZ13': '3HZ1',
        'HZ21': '1HZ2',
        'HZ22': '2HZ2',
        'HZ23': '3HZ2',

        'OXT': 'OXT',
        'H1' : '1H',
        'H2':  '2H',
        'H3':  '3H'
    }

    return amber_to_rosetta_atom_names.get(key)


#
def _map_atoms_pose_to_struct(key):

    # TODO: this should really be a generalized BCL substructure match at the residue level in C++ Rosetta
    # to create maps on-the-fly;
    # Rosetta --> Amber backbone atom map
    rosetta_to_amber_atom_names = {
        'N': 'N',
        'H': 'H',
        'CA': 'CA',
        'C': 'C',
        'O': 'O',
        '1HA': 'HA2',
        '2HA': 'HA3',

        '3HB': 'HB1',
        '1HB': 'HB2',
        '2HB': 'HB3',

        '3HG': 'HG1',
        '1HG': 'HG2',
        '2HG': 'HG3',

        '3HD': 'HD1',
        '1HD': 'HD2',
        '2HD': 'HD3',

        '3HE': 'HE1',
        '1HE': 'HE2',
        '2HE': 'HE3',

        # Branched sidechains


        'OXT': 'OXT'
        # '1H' : 'H1',
        # '2H': 'H2',
        # '3H': 'H3',
        # 'CN': 'CH3',
        # 'NM': 'N',
        # 'CO': 'C',
        # 'OP1': 'O',
        # 'CP2': 'CH3'
    }

    return rosetta_to_amber_atom_names.get(key)


#
def _get_rosetta_ace_cap_atoms():
    return ['CO', 'OP1', 'CP2', '1HP2', '2HP2', '3HP2']


#
def _get_rosetta_nme_cap_atoms():
    return ['NM', 'HM', 'CN', '1HN', '2HN', '3HN']


#
def _map_atoms_ace_pose_to_struct(key):
    rosetta_to_amber_atom_names = {
        'CO': 'C',
        'OP1': 'O',
        'CP2': 'CH3',
        '1HP2': 'H1',
        '2HP2': 'H2',
        '3HP2': 'H3'
    }

    return rosetta_to_amber_atom_names.get(key)


#
def _map_atoms_nme_pose_to_struct(key):
    rosetta_to_amber_atom_names = {
        'NM': 'N',
        'HM': 'H',
        'CN': 'CH3',
        '1HN': 'H1',
        '2HN': 'H2',
        '3HN': 'H3'
    }

    return rosetta_to_amber_atom_names.get(key)


#
def _map_resnames_three_to_one(key):
    # Define a lookup table relating the 3 and 1 letter codes
    # Is this really nowhere else? This should be somewhere accessible, I would think?
    # TODO: if not accessible somewhere, move to some common namespace functions file
    d = {
        'GLY' : 'G',
        'CYS' : 'C',
        'ASP' : 'D',
        'SER' : 'S',
        'GLN' : 'Q',
        'LYS' : 'K',
        'ILE' : 'I',
        'THR' : 'T',
        'PHE' : 'F',
        'ASN' : 'N',
        'HIS' : 'H',
        'LEU' : 'L',
        'ARG' : 'R',
        'TRP' : 'W',
        'ALA' : 'A',
        'VAL' : 'V',
        'GLU' : 'E',
        'TYR' : 'Y',
        'PRO' : 'P',
        'MET' : 'M',
        'HID' : 'H[HIS_D]',
        'ASH' : 'D[ASP_P1]',
        'CYM' : 'C[CYZ]',
        'CYX' : 'C[CYX]'
    }

    return d.get(key)


#
def _map_resnames_amber_to_rosetta(key):
    # Define a lookup table relating the 3 and 1 letter codes
    # Is this really nowhere else? This should be somewhere accessible, I would think?
    # TODO: if not accessible somewhere, move to some common namespace functions file
    d = {
        'ASH' : 'ASH',
        'CYM' : 'CYM',
        'CYX' : 'CYX',
        'HIE' : 'HIS',
        'HID' : 'HID'
    }

    return d.get(key)


#
def _set_atom_coords_struct_to_pose( pose, res, resid, struct, struct_resid, offset=0):
    struct_res = struct.residues
    for atnum in range(1, len(res.atoms()) + 1):
        matched_pose_atom = False
        atname = res.atom_name(atnum).strip()
        # ... and check against the atom names of the corresponding Structure residue on the inner loop
        for struct_atom in struct_res[struct_resid + offset].atoms:
            struct_resname = struct_res[struct_resid + offset].name
            struct_atname = struct_atom.name.strip()
            if struct_atname == atname:
                try:
                    struct_coords = \
                        Vec3D(struct.atoms[struct_atom.idx].xx, \
                              struct.atoms[struct_atom.idx].xy, \
                              struct.atoms[struct_atom.idx].xz)
                    # Set the Pose coordinates to the coordinates of the Structure atoms
                    pose.set_xyz(AtomID(atnum, resid), struct_coords)
                    matched_pose_atom = True
                    break
                except:
                    raise RosettaError('Could not set coordinates.')
        # Increment Pose atom
        atnum += int(1)
        # if not matched_pose_atom:
        #     print("WARNING: Atom " + \
        #           str(resid) + " " + str(res.name()) + " " + atname + \
        #           " has no atom name match from new Pose")


#
def _generate_seq_struct_to_pose( struct):
    pose_str = ""
    for resid in range(len(struct.residues)):
        if str(struct.residues[resid].name) == "ACE":
            pose_str += (str(_map_resnames_three_to_one(struct.residues[resid+1].name)) + "[" + struct.residues[
                resid+1].name + ":AcetylatedNtermProteinFull")
            if str(struct.residues[resid+2].name) == "NME":
                pose_str += ":MethylatedCtermProteinFull]"
            else:
                pose_str += "]"
        elif str(struct.residues[resid].name) == "NME" and str(struct.residues[resid-2].name) != "ACE":
            pose_str += (str(_map_resnames_three_to_one(struct.residues[resid-1].name)) + "[" + struct.residues[
                resid-1].name + ":MethylatedCtermProteinFull]")
        elif str(struct.residues[resid-1].name) == "ACE":
            continue
        elif _map_resnames_amber_to_rosetta(struct.residues[resid].name) is not None:
            three_letter_name = _map_resnames_amber_to_rosetta(struct.residues[resid].name)
            pose_str += str(_map_resnames_three_to_one(three_letter_name))
            print(str(three_letter_name) + " --> " + str(_map_resnames_three_to_one(three_letter_name)))
        elif _map_resnames_three_to_one(struct.residues[resid].name) is not None:
            pose_str += str(_map_resnames_three_to_one(struct.residues[resid].name))
    return pose_str


#
def _add_atoms_to_struct_from_pose(struct, pose, conf, res, resid, resname, atnum, offset=0, autodetect_bonds=False):
    chain = chr(res.chain() + ord('A') - 1)
    for atno, at in enumerate(res.atoms(), start=1):
        try:
            # Get atom info
            atinfo = res.atom_type(atno)
            atname = res.atom_name(atno).strip()

            # Debug
            # print("RESNAME - RESID - ATNAME: " + str(resname) + " - " + str(resid) + " - " + str(atname))

            # Rename ACE atoms
            if resname == "ACE":
                if atname not in _get_rosetta_ace_cap_atoms():
                    atnum += 1
                    continue
                elif _map_atoms_ace_pose_to_struct(atname) is not None:
                    atname = _map_atoms_ace_pose_to_struct(atname)

            # Rename NME atoms
            elif resname == "NME":
                if atname not in _get_rosetta_nme_cap_atoms():
                    atnum += 1
                    continue
                elif _map_atoms_nme_pose_to_struct(atname) is not None:
                    atname = _map_atoms_nme_pose_to_struct(atname)

            # Rename all other atoms
            elif _map_atoms_pose_to_struct(atname) is not None:
                atname = _map_atoms_pose_to_struct(atname)

            # If we are trying to do atoms in the base residue but come up with a cap atom then skip
            if resname != "ACE" and atname in _get_rosetta_ace_cap_atoms():
                atnum += 1
                continue
            elif resname != "NME" and atname in _get_rosetta_nme_cap_atoms():
                atnum += 1
                continue

            # If virtual atom handle with care
            if atinfo.is_virtual():
                atsym = 'EP'
            else:
                atsym = atinfo.element()
            rmin = atinfo.lj_radius()
            epsilon = atinfo.lj_wdepth()
            atomic_number = AtomicNum[atsym]
            mass = Mass[atsym]
        except KeyError:
            raise RosettaError('Could not recognize element: %s.'
                               % atsym)
        params = dict(atomic_number=atomic_number, name=atname,
                      charge=0.0, mass=mass, occupancy=0.0,
                      bfactor=0.0, altloc='', number=atnum,
                      rmin=rmin, epsilon=epsilon)
        if atinfo.is_virtual():
            atom = ExtraPoint(**params)
        else:
            atom = Atom(**params)

        # Set coordinates of new struct atom from pose
        atom.xx, atom.xy, atom.xz = (at.xyz()[0], at.xyz()[1], at.xyz()[2])

        struct.add_atom(atom, resname, resid+offset, chain, '')
        # struct.save("blah.add_atoms." + str(atnum)+".pdb")
        atnum += 1
        try:
            for nbr in conf.bonded_neighbor_all_res(AtomID(atno,resid)):
                # print("nbr.rsd(): " + str(nbr.rsd()) + ", resid: " + str(resid))
                # print("nbr.atomno(): " + str(nbr.atomno()) + ", atno: " + str(atno))
                # print("not autodetect_bonds: " + str(not autodetect_bonds))
                # print("atom: " + str(atom))
                if nbr.rsd() < resid or (nbr.rsd() == resid and nbr.atomno() < atno and not autodetect_bonds):
                    struct.bonds.append(Bond(struct.atoms[_n_prior(pose, nbr)],atom))
        except:
            raise RosettaError('Could not add bonds.')


#
def _parse_terminal_caps_pose_to_struct( struct, pose, resid, res_fullname_list):
    print("do stuff")


#
class RosettaPose(object):

    @staticmethod
    def load(pose):
        """
        Load a :class:`Pose` object and return a populated :class:`Structure`
        instance

        Parameters
        ----------
        pose : :class:`Pose`
            PyRosetta :class:`Pose` object to convert
        """
        if not Pose or not AtomID:
            raise ImportError('Could not load the PyRosetta module.')
        if not isinstance(pose, Pose):
            raise TypeError('Object is not a PyRosetta Pose object.')

        struct = Structure()

        atnum = 1
        conf = pose.conformation()
        for resid in range(1, pose.total_residue()+1):
            # Convenience
            res = pose.residue(resid)

            # To handle cap patches we need to parse the full residue name
            # print("Pose residue name: " + str(res.name()))
            res_fullname_list = res.name().split(":")
            # print("res fullname: " + str(res_fullname_list))
            # print( 'check: ' + str(res_fullname_list[0]) + " : " + str(res.name3().strip()))
            if res_fullname_list[0] == "HIS_D":
                resname = "HID"
            elif res_fullname_list[0] == "ASP_P1":
                resname = "ASH"
            elif res_fullname_list[0] == "CYZ":
                resname = "CYM"
            else:
                resname = res.name3().strip()

            # TODO: what follows is an initial parsing of terminal patches in Rosetta, but ultimately
            # this should probably be managed by a separate class. Will revisit this later hopefully
            # after learning a bit more about the intended class designs for objects that communicate
            # with ParmEd Structure objects from the ParmEd/OpenMM developer groups

            # Standard simple case where the res.name3 matches the full res.name
            # I also include the standard N- and C-term caps (for the peptide backbone zwitterion)
            # because those are auto detected easily and the atom types are easy to map
            # This is basically what the original implementation did
            if \
            (len(res_fullname_list) == int(1) and res_fullname_list[0] == res.name3()) or \
            (len(res_fullname_list) == int(1) and res_fullname_list[0] == "HIS_D") or \
            (len(res_fullname_list) == int(1) and res_fullname_list[0] == "ASP_P1") or \
            (len(res_fullname_list) == int(1) and res_fullname_list[0] == "CYZ") or \
            (len(res_fullname_list) == int(1) and res_fullname_list[0] == "CYX") or \
            (len(res_fullname_list) > int(1) and res_fullname_list[1] == "NtermProteinFull") or \
            (len(res_fullname_list) > int(1) and res_fullname_list[1] == "CtermProteinFull"):
                # Add atoms to struct based on pose atom type
                # print("RES_FULLNAME_LIST: " + str(res_fullname_list))
                # print("RESNAME: " + str(resname))
                _add_atoms_to_struct_from_pose( struct, pose, conf, res, resid, resname, atnum)

            # Now the case where we have both ACE and NME caps on a single residue
            # This is the classic dipeptide scenario
            elif \
            len(res_fullname_list) == int(3) and \
            res_fullname_list[0] == res.name3() and \
            ((res_fullname_list[1] == "AcetylatedNtermProteinFull" and res_fullname_list[2] == "MethylatedCtermProteinFull") or \
            (res_fullname_list[1] == "MethylatedCtermProteinFull" and res_fullname_list[2] == "AcetylatedNtermProteinFull")):
                _add_atoms_to_struct_from_pose( struct, pose, conf, res, resid, "ACE", atnum, offset=0, autodetect_bonds=True)
                _add_atoms_to_struct_from_pose(struct, pose, conf, res, resid, resname, atnum, offset=1, autodetect_bonds=True)
                _add_atoms_to_struct_from_pose(struct, pose, conf, res, resid, "NME", atnum, offset=2, autodetect_bonds=True)
                struct.assign_bonds()

        # Finalize
        struct.strip('@V1')  # Remove virtual atoms: CYX
        struct.strip('@NV')  # Remove virtual atoms: PRO?
        struct.strip('@VHG') # Remove virtual atoms: PRO?
        struct.unchange()
        return struct


    @staticmethod
    def struct_to_pose( struct):
        """
        Load a :class:`Structure` object and return a populated :class:`Pose`
        instance

        Parameters
        ----------
        struct : :class:`Structure`
            ParmEd :class:`Structure` object to convert
        reorder : ``bool``
            If true, reorder the resultant pose such that the sequence of atom names matches
            that of the input pose; otherwise, adopt the Rosetta pose atom name ordering
            (default=False)
        """
        if not Pose or not AtomID:
            raise ImportError('Could not load the PyRosetta module.')
        if not isinstance(struct, Structure):
            raise TypeError('Object is not a ParmEd Structure object.')

        # Initialize an empty pose
        pose = Pose()

        # I suspect that users going from MM to Rosetta will want full atom
        # Consider allowing this to be an option later and/or integrating non-standard all-atom types
        rts = rosetta.core.chemical.ChemicalManager.get_instance().residue_type_set('fa_standard')

        # The most irritating complication is that Rosetta treats NME and ACE caps just as additional
        # patched atoms on the base residue instead of as separate residues.
        # Check for ACE and NME caps; do not assume the entire sequence is one entire chain
        ace = []
        nme = []
        for res in range(len(struct.residues)):
            if str(struct.residues[res].name) == "ACE":
                ace.append(res + 1) # Rosetta indexing - before residue, so set to residue index by moving forward 1
            if str(struct.residues[res].name) == "NME":
                nme.append(res - 1) # Rosetta indexing - after residue, so set to residue index by moving backward 1

        # Build the sequence allowing autodetection of N- and C-terminal cap atoms
        pose_str = _generate_seq_struct_to_pose( struct)
        print(pose_str)
        rosetta.core.pose.make_pose_from_sequence(pose, pose_str, rts, True, True)

        # Map some atom names
        for atom in struct.atoms:
            mapped_name = _map_atoms_struct_to_pose(atom.name, atom.residue.name)
            if mapped_name is not None and atom.residue.name != "ACE" and atom.residue.name != "NME":
                atom.name = mapped_name
            else:
                mapped_name = _map_cap_atoms_struct_to_pose(atom.name, atom.residue.name)
                if mapped_name is not None:
                    atom.name = mapped_name

        # Now we need to set the coordinate positions of our atoms in our pose based on the original Structure
        # We should take care to preserve the atom ordering of our Structure
        struct_atoms = struct.atoms
        struct_res = struct.residues
        struct_resid = int( 0)
        # Consider refactoring to use a map of Structure Atom objects to Rosetta AtomID objects
        # Loop over all residues in Pose
        for resid in range(1, pose.total_residue() + 1):
            res = pose.residue(resid)
            resname = res.name3().strip()
            # Inefficient, but unlikely to be the bottleneck in an MM or design workflow
            # Iterate over the Pose atom names for current residue in outer loop...
            _set_atom_coords_struct_to_pose( pose, res, resid, struct, struct_resid)
            if resid in ace:
                _set_atom_coords_struct_to_pose( pose, res, resid, struct, struct_resid, 1)
            if resid in nme:
                _set_atom_coords_struct_to_pose( pose, res, resid, struct, struct_resid, -1)
            # Increment Pose residue
            struct_resid += int( 1)
        return pose