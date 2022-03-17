
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
import logging

multi_process_num = 3
write_step = 30

def get_logger(process_index):
    """get logger"""
    logger = logging.getLogger()  #获取logger
    logger.setLevel(logging.INFO)  # 设置记录等级时debug
    log_file = "{}.log".format(process_index)
    fh_all = logging.FileHandler(log_file, mode='w', encoding='utf-8')
    fh_all.setLevel(logging.INFO)
    fmt1 = logging.Formatter(fmt="%(asctime)s - %(levelname)-9s - %(filename)-8s : %(lineno)s line - %(message)s")
    fh_all.setFormatter(fmt1)
    logger.addHandler(fh_all)
    # logger.error("[error] input")
    return logger


def start_filter_multi_process(input_dir, file_lst, output_dir, process_index):
    # #####################################################################################
    logger = get_logger(process_index)
    # #####################################################################################
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
        atom_num = len(struct)
        print(atom_num)
        for isite in range(len(lse.coordination_environments)):
            type = lse.coordination_environments[isite][0]["ce_symbol"]
            print(type)
            if type != "T:4":
                return False
        return True
    for i, file in enumerate(file_lst):
        try:
            # if file != "Ba3(SnP2)2-mp-601867.cif":
            #     continue
            if i < 11460 // multi_process_num:
                continue
            if i % write_step == 0:
                logger.error(i)
            print("index = ", i, "file = ", file)
            address = input_dir + file
            if judge_tetrahedron_csm(address):
                print(file, "copy")
                shutil.copyfile(address, output_dir + file)
        except Exception as e:
            print(e)
            print(file, " error!!!")
        # break


if __name__ == "__main__":
    file_lst = os.listdir("../data/unit_cell/")
    file_lst.sort(key=lambda x:x[0])
    # file_lst = os.listdir("../data/tetrahedron_126335_csm/")

    file_lst_lst = [[] for _ in range(multi_process_num)]
    # 11460
    for k, file in enumerate(file_lst):
        if k < 60000:
            file_lst_lst[k % multi_process_num].append(file)
    for i in range(multi_process_num):
        print(len(file_lst_lst[i]))

    for i in range(multi_process_num):
        p = Process(target=start_filter_multi_process, args=("../data/unit_cell/", file_lst_lst[i],
                                               "../data/tetrahedron_126335_csm/", i),
                    name=f"work_{i}")
        p.start()
