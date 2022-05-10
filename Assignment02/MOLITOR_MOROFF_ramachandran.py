"""Dominik Molitor, Friederike Moroff, 12.05.2022, Assignment 01"""
import argparse
import math

from Bio.PDB import *
import matplotlib.pyplot as plt
import numpy as np

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
    # x = (v.x * w.x + v.y * w.y + v.z * w.z) / (math.sqrt(v.x ** 2 + v.y ** 2 + v.z ** 2) * math.sqrt(w.x ** 2 + w.y ** 2 + w.z ** 2))

    a = (v.x * w.x + v.y * w.y + v.z * w.z)
    b = math.sqrt(v.x ** 2 + v.y ** 2 + v.z ** 2)
    c = math.sqrt(w.x ** 2 + w.y ** 2 + w.z ** 2)
    x = a / (b * c)


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


def parse_file(id, file):
    parser = PDBParser()
    structure = parser.get_structure(id, file)

    return structure


def extract_coordinates(structure):
    all_atoms = []

    model = structure.get_models()
    models = list(model)

    for model in models:
        chain = model.get_chains()
        chains = list(chain)

        for chain in chains:
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


def compute_psi(atoms):

    for i in range(0, len(atoms)):
        if (atoms[i].name.strip() == "N" and atoms[i + 1].name.strip() == "CA" and atoms[i + 2].name.strip() == "C"):
            first_n = atoms[i]
            ca = atoms[i + 1]
            c = atoms[i + 2]
            for i in range(i + 1, len(atoms)):
                if (atoms[i].name.strip() == "N"):
                    second_n = atoms[i]

                    first_vector = get_normal_vector(first_n.vector, ca.vector, c.vector)
                    second_vector = get_normal_vector(ca.vector, c.vector, second_n.vector)

                    vector_angle(first_vector, second_vector)

    x = 0


def compute_phi(atoms):
    for i in range(0, len(atoms)):
        if (atoms[i].name.strip() == "C"):
            first_c = atoms[i]
            for i in range(i + 1, len(atoms)):
                if (atoms[i].name.strip() == "N" and atoms[i + 1].name.strip() == "CA" and atoms[i + 2].name.strip() == "C"):
                    n = atoms[i]
                    ca = atoms[i + 1]
                    second_c = atoms[i + 2]

                    first_vector = get_normal_vector(first_c.vector, n.vector, ca.vector)
                    second_vector = get_normal_vector(n.vector, ca.vector, second_c.vector)

                    vector_angle(first_vector, second_vector)


def plot_ramachandran (phi, psi):

    fig, ax = plt.subplots()

    ax.scatter(phi, psi)
    ax.set(xlim = (-180,180),
           ylim = (-180,180))
    plt.hlines(y = 0, xmin = -180, xmax = 180, colors= "grey")
    plt.vlines(x = 0, ymin = -180, ymax = 180, colors= "grey")
    plt.show()





if __name__ == '__main__':
    ## parsing input
    #parser = argparse.ArgumentParser()
    #parser.add_argument('-i', '--input_file', type=argparse.FileType('r'), nargs='+')
    #parser.add_argument('-o', '--output_file', type=str) # .pdf!!

    #args = parser.parse_args()

    #for f in args.input_file:
    #    print(f)

    #print(args.output_file)

    file_name = "/Users/friederike/Documents/Universität/Bioinformatik_Master/3_semester/structure_systems/assignments/SSBI/Assignment02/1igt.pdb"
    structure = parse_file("igt", file_name)
    all_atoms = extract_coordinates(structure)

    # compute_phi(all_atoms)
    compute_psi(all_atoms)

    p1 = Vector(1,2,3)
    p2 = Vector(4,5,6)

    angle = vector_angle(p1, p2)
    print(angle)
    # get_normal_vector(p1, p2, p3)


