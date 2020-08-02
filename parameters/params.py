# This script is to contain all the parameters within the simulation 

# from parameters.config import read_config
import yaml

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
Hilderbrand = {
       "PB":{"PS":  0.27062786*1e6},
       "PS":{"PB":  0.27062786*1e6,
             "PMMA":0.03006976*1e6},
       "PMMA":{"PS":0.03006976*1e6}
}

V_mono = {
       "PS":0.179*1e-27*6.02214086e23,
       "PB":0.111*1e-27*6.02214086e23,
       "PMMA":0.149*1e-27*6.02214086e23
}

SPECIES = params_dict["MATERIAL_CHOICE"]
chi_AB = Hilderbrand[SPECIES[0]][SPECIES[1]]*(V_mono[SPECIES[0]]*V_mono[SPECIES[1]])**0.5/params_dict["TEMPERATURE"]
N_A = params_dict["N_A"]
N_B = params_dict["N_B"]
# Insert other declarations 