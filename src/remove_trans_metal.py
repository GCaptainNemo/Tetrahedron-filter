from pymatgen.io.cif import CifParser
import os
import shutil


def is_transiant_metal(dir):
    parser = CifParser(dir)
    struct = parser.get_structures()[0]
    sites_lst = struct.sites
    element_lst = [site.specie for site in sites_lst]
    print(struct, struct[0])
    for i, element in enumerate(element_lst):
        if element.is_transition_metal and element.Z not in [29, 30, 47, 48, 79, 80]:
            return True
    return False
   

if __name__ == "__main__":
    dir = "../data/tetrahedron_126335_csm_r_repeat_seperate_1237/"
    for i in range(1, 6):
        file_dir = dir + "{}/".format(i)
        new_file_dir = dir + "{}_remove_transiant_metal/".format(i)
        if not os.path.exists(new_file_dir):
            os.mkdir(new_file_dir)
        file_lst = os.listdir(file_dir)
        for j, file in enumerate(file_lst):
            if not is_transiant_metal(file_dir + file):
                shutil.copy(file_dir + file, new_file_dir + file)