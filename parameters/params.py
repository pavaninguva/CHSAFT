# This script is to contain all the parameters within the simulation 

# from parameters.config import read_config
import yaml
from Thermo.Thermo import GibbsMixingFH

with open("params.yml", 'r') as stream:
    params_dict = yaml.safe_load(stream)
    print (params_dict)

# Cahn Hilliard Model type: 
MOBILITY_MODEL = params_dict["MOBILITY_MODEL"]

# Compositon parameters
A_RAW = params_dict["A_RAW"]
NOISE_MAGNITUDE = params_dict["NOISE_MAGNITUDE"]

# Material parameters (either PS_PMMA / PS_PB)
MATERIAL_CHOICE = params_dict["MATERIAL_CHOICE"]
SIZE_DISPARITY = params_dict["SIZE_DISPARITY"]


# Homogenous Free energy function (PC_SAFT, FH, etc)
GIBBS = params_dict["GIBBS"]
TEMPERATURE = params_dict["TEMPERATURE"]

# Numerics
DT = params_dict["DT"]
TIME_MAX = params_dict["TIME_MAX"]
N_CELLS = params_dict["N_CELLS"]
DOMAIN_LENGTH = params_dict["DOMAIN_LENGTH"]
theta_ch = params_dict["theta_ch"]
MESH_TYPE = params_dict["MESH_TYPE"]
TIME_STRIDE = params_dict["TIME_STRIDE"]
FINITE_ELEMENT = params_dict["FINITE_ELEMENT"]
FINITE_ELEMENT_ORDER = params_dict["FINITE_ELEMENT_ORDER"]
SOLVER_CONFIG = params_dict["SOLVER_CONFIG"]

# Logic for material interaction parameters
A_SYM = []
if MATERIAL_CHOICE == "PS_PMMA":
    N_A = params_dict["N_A_PS_PMMA"]
    N_B = params_dict["N_B_PS_PMMA"]
    SPECIES = ["PS","PMMA"]
    r = GibbsMixingFH(SPECIES,[N_A,N_B],[TEMPERATURE])
    chi_AB = r.chi_AB()
    # Insert other declarations  
elif MATERIAL_CHOICE == "PS_PB":
    N_A = params_dict["N_A_PS_PB"]
    N_B = params_dict["N_B_PS_PB"]
    SPECIES = ["PS","PB"]
    r = GibbsMixingFH(SPECIES,[N_A,N_B],[TEMPERATURE])
    chi_AB = r.chi_AB()
    # Insert other declarations 
elif MATERIAL_CHOICE == "Vashistha94":
    chi_AB = params_dict["chi_AB_Vashistha"]
    N_A = params_dict["N_A_Vashistha"]
    N_B = params_dict["N_B_Vashistha"]
    SPECIES = ["Vashistha94"]
    A_SYM = params_dict["A_SYM_Vashistha"]
elif MATERIAL_CHOICE == "He97":
    chi_AB = params_dict["chi_AB_He"]
    N_A = params_dict["N_A_He"]
    N_B = params_dict["N_B_He"]
    SPECIES = ["He97"]
    A_SYM = params_dict["A_SYM_He"]
# Insert other declarations 

