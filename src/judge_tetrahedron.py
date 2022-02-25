#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author： 11360
# datetime： 2021/6/2 12:20 

import ase
from ase.io import read
import os
import numpy as np
import shutil
import pickle

from pymatgen.core.periodic_table import Element

class ElementData:
    # given formula return atomic number, e.g.,  element_dict["H"] = 1
    element_dict = ase.data.atomic_numbers

    # formula name dict
    formula_name_lst = list(ase.data.atomic_numbers.keys())

def is_angles(vector1, vector2, precision=12):
    """
    判断两向量夹角是不是四面体夹角 109°
    """
    cos_angle = vector1 @ vector2 / np.linalg.norm(vector1) / np.linalg.norm(vector2)

    angle = np.arccos(cos_angle) * 180 / np.pi
    print("angle = ", angle)
    # if np.abs(angle - 109) < precision:
    if 95 < angle < 150:
        return True
    else:
        return False


def extend_atom(atoms_obj):
    """
    对晶胞进行扩胞处理(3x3x3)
    """
    pos = atoms_obj.get_positions()
    cell = atoms_obj.get_cell()
    new_atoms_obj = atoms_obj.copy()
    # 12
    vec_a = cell[0]
    vec_b = cell[1]
    vec_c = cell[2]
    for k in range(-1, 2, 1):
        for i in range(2):
            pos_a = pos + cell[i] + vec_c * k
            new_atoms_obj.set_positions(newpositions=pos_a)
            atoms_obj.extend(new_atoms_obj)

            pos_a_n = pos - cell[i] + vec_c * k
            new_atoms_obj.set_positions(newpositions=pos_a_n)
            atoms_obj.extend(new_atoms_obj)
    # 2
    for i in range(-1, 2, 2):
        pos_a = pos + vec_c * i
        new_atoms_obj.set_positions(newpositions=pos_a)
        atoms_obj.extend(new_atoms_obj)

    # 12
    for k in range(-1, 2, 1):
        for i in range(2):
            for j in range(2):
                new_pos = pos + (-1) ** i * vec_a + (-1) ** j * vec_b + k * vec_c
                new_atoms_obj.set_positions(newpositions=new_pos)
                atoms_obj.extend(new_atoms_obj)
    return atoms_obj


def judge_tetrahedron_bond(atoms_obj, origin_num, tolerance=0.55):
    """
    对扩胞后判断是否为四面体,使用成键经验公式
    """
    # ############################################################
    # judge whether bonding

    atomic_numbers = atoms_obj.get_atomic_numbers()
    dist_matrix = atoms_obj.get_all_distances(vector=False)
    dist_matrix += np.diag([np.inf for i in range(dist_matrix.shape[0])])  # 不考虑对角元素
    # print("dist_matrix = ", dist_matrix)
    ion_radii_vec = [Element(ElementData.formula_name_lst[atomic_num]).average_ionic_radius
                     for atomic_num in atomic_numbers]
    ion_radii_vec = np.reshape(np.array(ion_radii_vec), [-1, 1])
    ion_radii_mat = (ion_radii_vec + np.transpose(ion_radii_vec)) * (1. + tolerance)
    diff_matrix = dist_matrix - ion_radii_mat
    # ############################################################
    for atom in range(origin_num):
        linshi = diff_matrix[atom, :]
        bond_index = np.where(linshi < 0)[0]
        num = bond_index.shape[0]
        if num != 4:
            # a = [138, 8, 9, 139]
            # for i in a:
            #     print(linshi[i])
            print("atom = {}, num = {} != 4".format(atom, num))
            return False
    return True


def judge_tetrahedron_4nn(atoms_obj, origin_num):
    """
    扩胞后,仅用4个邻居，几何结构进行判断，未使用成键经验公式 128/22377
    """
    # ############################################################
    dist_matrix = atoms_obj.get_all_distances(vector=False)[:origin_num, :]
    dist_vec = atoms_obj.get_all_distances(vector=True)[:origin_num, :, :]
    print(dist_vec.shape)
    # ############################################################
    # use bonding angle, length
    # ############################################################
    for atom in range(origin_num):
        print("***********")
        print("atom:", atom)
        linshi = dist_matrix[atom, :]
        index = linshi.argsort()[:5]
        print("dis = ", linshi[index])
        _ = np.argmin(linshi[index])
        index = np.delete(index, _)
        dis = linshi[index]
        print("dis = ", dis)
        print("index = ", index)
        if np.max(dis).item() - np.min(dis).item() > 0.4:
            print("bond length unequal")
            return False
        linshi_vec = dist_vec[atom, :]
        for a1 in range(4):
            for a3 in range(a1 + 1, 4):
                vec1 = linshi_vec[index[a1]]
                print("vec1 = ", vec1)
                vec2 = linshi_vec[index[a3]]
                print("vec2 = ", vec2)
                if not is_angles(vec1, vec2):
                    print("bond angle != 109")
                    return False
    print("return True")
    return True


def start_filter(input_dir="../data/conventional_cell/",
                 judge_function=judge_tetrahedron_4nn,
                 output_dir="../data/tetrahedron_124657/"):
    if not os.path.exists(output_dir):
        raise ValueError("[Error] output dir does not exist")
    too_big_file = []
    file_lst = os.listdir(input_dir)
    for i, file in enumerate(file_lst):
        # if file != "Cu4 B4 S8-mp-12954.cif":
        #     continue
        # if i < 690:
        #     continue
        try:
            print("index = ", i)
            address = input_dir + file
            atoms_obj = ase.io.read(address)
            origin_num = len(atoms_obj)
            if origin_num > 200:
                too_big_file.append(file)
                continue

            atoms_obj = extend_atom(atoms_obj)
            print(" ------------------- ")
            print(file)
            if judge_function(atoms_obj, origin_num):
                print("copy")
                shutil.copyfile(address, output_dir + file)
        except Exception as e:
            print(e)
            print(file, " error!!!")
        # break
    print(too_big_file)

    with open("too_big_file.pkl", "wb") as f:
        pickle.dump(too_big_file, f)


if __name__ == "__main__":
    start_filter(input_dir="../data/conventional_cell/", judge_function=judge_tetrahedron_bond,
                 output_dir="../data/tetrahedron_124657_bond/")


