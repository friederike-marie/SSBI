"""Dominik Molitor, Friederike Moroff, 12.05.2022, Assignment 01"""
import argparse
import math

import numpy as np
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


# computes the angle between two vectors
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


# computes the normal of a pane given three point of it
def get_normal(u, v, w):
    # build the two direction vectors defining a plane
    d1 = Vector(v.x - u.x, v.y - u.y, v.z - u.z)
    d2 = Vector(w.x - u.x, w.y - u.y, w.z - u.z)

    # compute the normal from the two direction vectors
    a = d1.y * d2.z - d1.z * d2.y
    b = d1.z * d2.x - d1.x * d2.z
    c = d1.x * d2.y - d1.y * d2.x

    normal = Vector(a, b, c)

    return normal


# checks if a point is on the same side of the pane as it's normal
def clockwise(n, s, p):
    a = (n.x * p.x + n.y * p.y + n.z * p.z) - (s.x * n.x + s.y * n.y + s.z * n.z)
    b = math.sqrt(n.x ** 2 + n.y ** 2 + n.z ** 2)
    x = a / b

    return x >= 0


# parses the file and returns the structure
def parse_file(id, file):
    parser = PDBParser()
    structure = parser.get_structure(id, file)

    return structure


# extracts all atom coordinates from a given chain of the protein
def extract_coordinates(chain):
    all_atoms = []

    # gets all residues of chain
    residue = chain.get_residues()
    residues = list(residue)

    # gets all atoms of all residues
    for residue in residues:
        atom = residue.get_atoms()
        atoms = list(atom)

        # appends all atoms of the residue to the list of atoms
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

                    # builds a normal for the pane describing N, Ca, C
                    first_vector = get_normal(first_n.vector, ca.vector, c.vector)
                    # builds a normal for the pane describing Ca, C, N'
                    second_vector = get_normal(ca.vector, c.vector, second_n.vector)

                    # computes the angle between both panes
                    angle = vector_angle(first_vector, second_vector)
                    if not clockwise(first_vector, first_n.vector, second_n.vector):
                        angle = -angle
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

                    print(first_c.name)
                    print(str(first_c.vector.x))
                    print(str(first_c.vector.y))
                    print(str(first_c.vector.z))

                    print(n.name)
                    print(str(n.vector.x))
                    print(str(n.vector.y))
                    print(str(n.vector.z))

                    print(ca.name)
                    print(str(ca.vector.x))
                    print(str(ca.vector.y))
                    print(str(ca.vector.z))

                    print(second_c.name)
                    print(str(second_c.vector.x))
                    print(str(second_c.vector.y))
                    print(str(second_c.vector.z))

                    # builds a normal for the pane describing C, N, Ca
                    first_vector = get_normal(first_c.vector, n.vector, ca.vector)
                    # builds a normal for the pane describing N, Ca, C
                    second_vector = get_normal(n.vector, ca.vector, second_c.vector)

                    # computes the angle between both panes
                    angle = vector_angle(first_vector, second_vector)
                    if not clockwise(first_vector, first_c.vector, second_c.vector):
                        angle = -angle
                    all_phi.append(angle)

    return all_phi


# computes all torsion angles of a given structure
def compute_all_angles(structure):
    all_phi = []
    all_psi = []

    model = structure.get_models()
    models = list(model)

    # gets all chains of a model
    for model in models:
        chain = model.get_chains()
        chains = list(chain)

        # extracs for all chains all atoms an computes the phi and psi angle
        for chain in chains:
            atoms = extract_coordinates(chain)
            phi = compute_phi(atoms)
            psi = compute_psi(atoms)
            if phi:
                all_phi.extend(phi)
                all_psi.extend(psi)

    plt1 = plot_ramachandran(all_phi[0:500], all_psi[0:500])
    plt2 = plot_ramachandran(all_phi[0:4000], all_psi[0:4000])
    plt3 = plot_ramachandran(all_phi[0:8000], all_psi[0:8000])
    plt4 = plot_ramachandran(all_phi, all_psi)

    # output_file = PdfPages('1mbn.pdf')
    # output_file.savefig(plt1)
    # output_file.savefig(plt2)
    # output_file.savefig(plt3)
    # output_file.savefig(plt4)
    # output_file.close()


# plots list of given phi and psi angles
def plot_ramachandran(phi, psi):
    heatmap, xedges, yedges = np.histogram2d(phi, psi, bins = 50)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    plt.clf()
    plt.imshow(heatmap.T, extent=extent, origin='lower')
    plt.show()

    #fig, ax = plt.subplots()

    #ax.scatter(phi, psi, s=0.5, color="black")
    #ax.set(xlim=(-180, 180),
     #      ylim=(-180, 180))
    #plt.hlines(y=0, xmin=-180, xmax=180, colors="grey")
    #plt.vlines(x=0, ymin=-180, ymax=180, colors="grey")
    #plt.xlabel("phi")
    #plt.ylabel("psi")

    return heatmap


# checks if the argument for the output file name ends with .pdf
def pdf_validator(astring):
    if not isinstance(astring, str):
        raise ValueError
    if not astring.endswith(".pdf"):
        raise ValueError
    return astring


if __name__ == '__main__':
    # parsing input
    # parser = argparse.ArgumentParser()
    # parser.add_argument('-i', '--INPUT_File', type=argparse.FileType('r'), nargs='+')
    # parser.add_argument('-o', '--OUTPUT_FILE', type=pdf_validator)
    # args = parser.parse_args()

    # for f in args.INPUT_FILE:
    # for all files!

    file_name = "/Users/friederike/Documents/Universität/Bioinformatik_Master/3_semester/structure_systems/assignments/SSBI/Assignment02/2hik.pdb"
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
