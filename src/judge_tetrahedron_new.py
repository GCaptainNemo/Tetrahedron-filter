
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import \
    SimplestChemenvStrategy, MultiWeightsChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.structure_environments \
import LightStructureEnvironments
from pymatgen.ext.matproj import MPRester
from pymatgen.io.cif import CifParser
import os
import shutil
from multiprocessing import Process

multi_process_num = 3


def start_filter(input_dir, file_lst, output_dir):
    lgf = LocalGeometryFinder()
    lgf.setup_parameters(centering_type='centroid', include_central_site_in_centroid=True)
    def judge_tetrahedron_csm(unit_cell_cif_address):
        parser = CifParser(unit_cell_cif_address)
        struct = parser.get_structures()[0]
        lgf.setup_structure(structure=struct)

        # se
        se = lgf.compute_structure_environments(maximum_distance_factor=1.41, only_cations=False)
        # strategy
        strategy = SimplestChemenvStrategy(distance_cutoff=1.4, angle_cutoff=0.3)
        lse = LightStructureEnvironments.from_structure_environments(strategy=strategy,
                                                                     structure_environments=se)
        print(len(struct))
        atom_num = len(struct)
        for isite in range(atom_num):
            type = lse.coordination_environments[isite][0]["ce_symbol"][0]
            print(type)
            if type != "T":
                return False
        return True
    for i, file in enumerate(file_lst):
        try:
            if i < 4946 // multi_process_num:
                continue
            print("index = ", i, "file = ", file)
            address = input_dir + file
            if judge_tetrahedron_csm(address):
                print("copy")
                shutil.copyfile(address, output_dir + file)
        except Exception as e:
            print(e)
            print(file, " error!!!")


if __name__ == "__main__":
    file_lst = os.listdir("../data/unit_cell/")
    file_lst_lst = [[] for _ in range(multi_process_num)]
    for k, file in enumerate(file_lst):
        file_lst_lst[k % 3].append(file)
    for i in range(multi_process_num):
        print(len(file_lst_lst[i]))

    for i in range(multi_process_num):
        p = Process(target=start_filter, args=("../data/unit_cell/", file_lst_lst[i],
                                               "../data/tetrahedron_126335_csm/"),
                    name=f"work_{i}")
        p.start()
