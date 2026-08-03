import json
import shutil
import os
from tqdm import tqdm
import sys

list_num = sys.argv[1]

with open(f'list', 'r') as file:
    file_list = [line.strip() for line in file.readlines()]

with open('t_range', 'r') as file:
    t_range = [line.strip() for line in file.readlines()]

with open(f'p_range_{list_num}', 'r') as file:
    p_range = [line.strip() for line in file.readlines()]

with open('simulation.json','r') as file:
    input_template = json.load(file)

# sc_info = "mof_list_unitcell.dat"
# sc_dict = {}
# with open(sc_info, "r") as file:
#     for line in file:
#         parts = line.strip().split()
#         structure_name = parts[0]
#         unitcell_numbers = list(map(int, parts[1:]))
#         sc_dict[structure_name] = unitcell_numbers

if not os.path.exists('calculation'): os.makedirs('calculation')

for struc in tqdm(file_list):
    struc_dir = os.path.join('calculation',struc)
    if not os.path.exists(struc_dir): os.makedirs(struc_dir)
    for T in t_range:
        T_dir = os.path.join(struc_dir,f'{T}K')
        if not os.path.exists(T_dir): os.makedirs(T_dir)
        for P in p_range:
            P_dir = os.path.join(T_dir,f'{P}Pa')
            if not os.path.exists(P_dir): os.makedirs(P_dir)

            input_tmp = input_template
            input_tmp["Systems"][0]["Name"] = struc
            # input_tmp["Systems"][0]["NumberOfUnitCells"] = sc_dict[struc]
            input_tmp["Systems"][0]["ExternalTemperature"] = float(T)
            input_tmp["Systems"][0]["ExternalPressure"] = float(P)

            with open(os.path.join(P_dir,'simulation.json'),'w') as file:
                json.dump(input_tmp,file,indent=4)

            shutil.copy('force_field.json', P_dir)
            shutil.copy('co2.json', P_dir)
            shutil.copy('raspa3_run.slurm', P_dir)
            shutil.copy('mbx.json', P_dir)

            shutil.copy(os.path.join('../structures',f'{struc}.cif'), P_dir)

