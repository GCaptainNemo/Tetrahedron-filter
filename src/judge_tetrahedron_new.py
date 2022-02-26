
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import \
    SimplestChemenvStrategy, MultiWeightsChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.structure_environments \
import LightStructureEnvironments
from pymatgen.ext.matproj import MPRester
from pymatgen.io.cif import CifParser
import os
import shutil


def judge_tetrahedron_csm(unit_cell_cif_address):
    parser = CifParser(unit_cell_cif_address)
    struct = parser.get_structures()[0]
    lgf = LocalGeometryFinder()
    lgf.setup_parameters(centering_type='centroid', include_central_site_in_centroid=True)
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


def start_filter(input_dir="../data/unit_cell/",
                 judge_function=judge_tetrahedron_csm,
                 output_dir="../data/tetrahedron_126335_csm/"):
    if not os.path.exists(output_dir):
        raise ValueError("[Error] output dir does not exist")
    file_lst = os.listdir(input_dir)
    for i, file in enumerate(file_lst):
        try:
            print("index = ", i)
            address = input_dir + file
            if judge_function(address):
                print("copy")
                shutil.copyfile(address, output_dir + file)
        except Exception as e:
            print(e)
            print(file, " error!!!")


if __name__ == "__main__":
    start_filter(input_dir="../data/tetrahedron_124657_union_unit_cell/", judge_function=judge_tetrahedron_csm,
                 output_dir="../data/tetrahedron_126335_csm/")
