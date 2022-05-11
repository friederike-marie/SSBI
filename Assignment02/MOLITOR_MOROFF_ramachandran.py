"""Dominik Molitor, Friederike Moroff, 12.05.2022, Assignment 01"""
import argparse
import math

from Bio.PDB import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


class Vector:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z


class Atom:
    def __init__(self, name, vector):
        self.name = name
        self.vector = vector


def vector_angle(v, w):
    a = (v.x * w.x + v.y * w.y + v.z * w.z)
    b = math.sqrt(v.x ** 2 + v.y ** 2 + v.z ** 2)
    c = math.sqrt(w.x ** 2 + w.y ** 2 + w.z ** 2)
    x = a / (b * c)

    if x > 1:
        x = 1
    if x < -1:
        x = -1

    angle_degree = math.degrees(math.acos(x))

    return angle_degree


def get_normal_vector(u, v, w):
    # build the two direction vectors defining a plane
    d1 = Vector(v.x - u.x, v.y - u.y, v.z - u.z)
    d2 = Vector(w.x - u.x, w.y - u.y, w.z - u.z)

    # compute the normal vector from the two direction vectors
    a = d1.y * d2.z - d1.z * d2.y
    b = d1.z * d2.x - d1.x * d2.z
    c = d1.x * d2.y - d1.y * d2.x

    normal_vector = Vector(a, b, c)

    return normal_vector


def clockwise(n, s, p):
    a = (n.x * p.x + n.y * p.y + n.z * p.z) - (s.x * n.x + s.y * n.y + s.z * n.z)
    b = math.sqrt(n.x ** 2 + n.y ** 2 + n.z ** 2)
    x = a / b

    return x >= 0


def parse_file(id, file):
    parser = PDBParser()
    structure = parser.get_structure(id, file)

    return structure


def extract_coordinates(chain):
    all_atoms = []

    residue = chain.get_residues()
    residues = list(residue)

    for residue in residues:
        atom = residue.get_atoms()
        atoms = list(atom)

        for i in range(0, len(atoms)):
            atom = atoms[i]
            new_atom = Atom(atom.fullname, Vector(atom.coord[0], atom.coord[1], atom.coord[2]))
            all_atoms.append(new_atom)

    return all_atoms


# computes the psi angle by computing the angle of the normal vector from two panes
def compute_psi(atoms):
    all_psi = []
    for i in range(0, len(atoms)):
        if atoms[i].name.strip() == "N" and atoms[i + 1].name.strip() == "CA" and atoms[i + 2].name.strip() == "C":
            first_n = atoms[i]
            ca = atoms[i + 1]
            c = atoms[i + 2]
            for i in range(i + 1, len(atoms)):
                if atoms[i].name.strip() == "N":
                    second_n = atoms[i]

                    first_vector = get_normal_vector(first_n.vector, ca.vector, c.vector)
                    second_vector = get_normal_vector(ca.vector, c.vector, second_n.vector)

                    vector_angle(first_vector, second_vector)

                    angle = vector_angle(first_vector, second_vector)
                    if not clockwise(first_vector, first_n.vector, second_n.vector):
                        angle = -angle
                    print(angle)
                    all_psi.append(angle)

    return all_psi


# computes the phi angle by computing the angle of the normal vector from two panes
def compute_phi(atoms):
    all_phi = []
    for i in range(0, len(atoms)):
        if atoms[i].name.strip() == "C":
            first_c = atoms[i]
            for i in range(i + 1, len(atoms)):
                if atoms[i].name.strip() == "N" and atoms[i + 1].name.strip() == "CA" and atoms[i + 2].name.strip() == "C":
                    n = atoms[i]
                    ca = atoms[i + 1]
                    second_c = atoms[i + 2]

                    first_vector = get_normal_vector(first_c.vector, n.vector, ca.vector)
                    second_vector = get_normal_vector(n.vector, ca.vector, second_c.vector)

                    angle = vector_angle(first_vector, second_vector)
                    if not clockwise(first_vector, first_c.vector, second_c.vector):
                        angle = -angle
                    print(angle)
                    all_phi.append(angle)

    return all_phi


def compute_all_angles(structure):
    all_phi = []
    all_psi = []

    model = structure.get_models()
    models = list(model)

    for model in models:
        chain = model.get_chains()
        chains = list(chain)

        for chain in chains:
            atoms = extract_coordinates(chain)
            phi = compute_phi(atoms)
            psi = compute_psi(atoms)
            if phi:
                all_phi.extend(phi)
                all_psi.extend(psi)
    print(len(all_phi))
    print(len(all_psi))

    plot_ramachandran(all_phi[0:500], all_psi[0:500])
    plot_ramachandran(all_phi[0:1000], all_psi[0:1000])
    plot_ramachandran(all_phi[0:2000], all_psi[0:2000])
    plot_ramachandran(all_phi[0:4000], all_psi[0:4000])
    plot_ramachandran(all_phi[0:8000], all_psi[0:8000])
    plot_ramachandran(all_phi, all_psi)


# plots list of given phi and psi angles
def plot_ramachandran(phi, psi):
    fig, ax = plt.subplots()

    ax.scatter(phi, psi, s=0.5, color="black")
    ax.set(xlim=(-180, 180),
           ylim=(-180, 180))
    plt.hlines(y=0, xmin=-180, xmax=180, colors="grey")
    plt.vlines(x=0, ymin=-180, ymax=180, colors="grey")
    plt.xlabel("phi")
    plt.ylabel("psi")
    plt.show()

    return fig


def pdf_validator(astring):
    if not isinstance(astring, str):
        raise ValueError
    if not astring.endswith(".pdf"):
        raise ValueError
    return astring


if __name__ == '__main__':
    # parsing input
    # parser = argparse.ArgumentParser()
    # parser.add_argument('-i', '--input_file', type=argparse.FileType('r'), nargs='+')
    # parser.add_argument('-o', '--output_file', type=pdf_validator)
    # args = parser.parse_args()

    # for f in args.input_file:
    # for all files!

    file_name = "/Users/friederike/Documents/Universität/Bioinformatik_Master/3_semester/structure_systems/assignments/SSBI/Assignment02/1igt.pdb"
    structure = parse_file("igt", file_name)

    compute_all_angles(structure)


    # output file -------------------------
    # p3 = [1,2,3]
    # p4 = [4,5,6]

    # plt1 = plot_ramachandran(p3, p4)
    # plt2 = plot_ramachandran(p3, p4)

    # output_file = PdfPages('output.pdf')
    # output_file.savefig(plt1)
    # output_file.savefig(plt2)
    # output_file.close()
