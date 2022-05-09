"""Dominik Molitor, Friederike Moroff, 12.05.2022, Assignment 01"""
import argparse
import math

from Bio.PDB import *
from numpy import arccos


class Vector:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z


def vector_angle(v, w):
    x = (v.x * w.x + v.y * w.y + v.z * w.z) / (math.sqrt(v.x ** 2 + v.y ** 2 + v.z ** 2) + math.sqrt(w.x ** 2 + w.y ** 2 + w.z ** 2))
    angle = math.degrees(math.acos(x))

    return angle


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
    atoms = structure.get_atoms()


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
    parse_file("igt", file_name)

    p1 = Vector(1, 0, 0)
    p2 = Vector(0, 1, 0)
    p3 = Vector(0, 0, 1)


    angle = vector_angle(p1, p2)

    get_normal_vector(p1, p2, p3)

