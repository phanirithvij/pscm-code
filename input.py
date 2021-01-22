"""
Read a pdb file and access frames with molecules, atoms
"""

import sys
from typing import List, Tuple


WARN_INVALID_MOLECULE = """
[WARNING]: Invalid Molecule
Try 1-based indexing not 0 based for molecules from frame
"""


class AtomBond():
    """An AtomBond, has an
        ```
        id idx: int,
        bond_name: str,
        postiion x, y, z float,
        molecule: Molecule
        ```
    """

    def __init__(self,
                 name: str,
                 idx: int,
                 x: float, y: float, z: float):
        self.id = idx
        self.molecule: Molecule = None
        self.x = x
        self.y = y
        self.z = z
        self.name = name

    def __repr__(self):
        return str(self.__dict__)

    def __str__(self):
        return f"<{self.__class__.__name__}@{self.id}>: {self.__dict__}"


class Molecule():
    """A Molecule, has `id idx: int, atoms: List[atoms]`"""

    def __init__(self, idx: int, atoms: List[AtomBond] = []):
        self.id = idx
        self.atoms = atoms

    def __repr__(self):
        return str(self.__dict__)

    def __str__(self):
        if self.id == 0:
            return str(WARN_INVALID_MOLECULE)
        return f"<{self.__class__.__name__}@{self.id}>: {self.__dict__}"


class Frame():
    """A frame from the pdb it has `list[atoms], list[molecules]`"""

    def __init__(self, num: int, atom_strs: List[str] = []):
        self.atoms: List[AtomBond] = []
        self.num = num
        # have an invalid molecule at 0 index
        # because input file has 1 based indexing
        self.molecules: List[Molecule] = [Molecule(0)]
        # ATOM   1296  H2  TIP3W 504     -10.719   2.716   2.401  0.00  0.00      W
        i = 0
        for a in atom_strs:
            i += 1
            typ, num, name, _, mol_id, x, y, z, _, _, _ = [
                x for x in a.split(' ') if x != '']
            mol_id = int(mol_id)
            num = int(num)
            x, y, z = float(x), float(y), float(z)
            if typ == 'ATOM':
                at = AtomBond(name, num, x, y, z)
                self.atoms.append(at)
                try:
                    at.molecule = self.molecules[mol_id]
                    self.molecules[mol_id].atoms.append(at)
                except IndexError:
                    # create a new molecule at mol_id
                    mol = Molecule(mol_id, [])
                    at.molecule = mol
                    mol.atoms.append(at)
                    self.molecules.insert(mol_id+1, mol)
            else:
                print("A non ATOM row found", a)
                exit(1)

    # Use repr if whole json is needed
    # Use str for human readable string
    # https://stackoverflow.com/a/1438297/8608146

    def __str__(self):
        dictx = {}
        for key in self.__dict__:
            if key == 'atoms':
                dictx['atoms'] = f"[{{...{len(self.atoms)} atoms...}}]"
            elif key == 'molecules':
                dictx['molecules'] = f"[{{...{len(self.molecules)} molecules...}}]"
            else:
                dictx[key] = self.__dict__[key]
        return f"<{self.__class__.__name__}@{self.num}>: {dictx}"

    def __repr__(self):
        return str(self.__dict__)


def read_pdb(filename: str) -> Tuple[bool, List[Frame]]:
    if not valid_pdb(filename):
        print("Invalid pdb file", filename)
        return False, []
    frames: List[Frame] = []
    atoms = []
    frame_no = 0
    with open(filename, 'r') as pdb:
        for i, line in enumerate(pdb.readlines()):
            if i == 0:
                # skip first line
                continue
            line = line.strip()
            if line == "END":
                f = Frame(frame_no, atoms)
                frames.append(f)
                atoms = []
                frame_no += 1
            else:
                atoms.append(line)
    return True, frames


def valid_pdb(filename: str):
    with open(filename, 'r') as pdb:
        for line in pdb.readlines():
            parts = [x for x in line.strip().split(' ') if x != '']
            # read the very first ATOM line and check if a valid line
            if parts[0] == 'ATOM':
                if len(parts) == 11:
                    return True
                break
        return False


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python this.py <filename.pdb>")
        exit(1)
    valid, frames = read_pdb(sys.argv[1])
    if valid:
        print(frames[0])
        print()
        print(frames[0].molecules[1])
        print()
        print(frames[0].atoms[0])
