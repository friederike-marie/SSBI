"""Dominik Molitor, Friederike Moroff, 12.05.2022, Assignment 01"""
import argparse
import math

from Bio.PDB import *
from numpy import arccos


class Point:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z


def vector_angle(v, w):
    x = (v.x * w.x + v.y * w.y + v.z * w.z) / (math.sqrt(v.x ** 2 + v.y ** 2 + v.z ** 2) + math.sqrt(w.x ** 2 + w.y ** 2 + w.z ** 2))
    angle = math.degrees(math.acos(x))

    return angle


def get_normal_vector(u, v, w):





def parse_file(id, file):
    parser = PDBParser()
    structure = parser.get_structure(id, file)
    print(structure.get_atoms())


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
    # parse_file("igt", file_name)

    point1 = Point(1, 0, 0)
    point2 = Point(0, 1, 0)
    point3 = Point(0, 0, 1)



    print(math.degrees(math.acos(0)))

    angle = vector_angle(point1, point2)
    print(angle)

