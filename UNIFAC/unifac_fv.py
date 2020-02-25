
from numpy import log as ln

"""
This section is for capturing the physical properties / group parameters of the polymers

We define the extent of polymerisation, number of groups etc
"""

# Update the input polymer properties as needed
Polymer_length = {
    "PMMA": 490,
    "PS":490,
    "PB":490
}

# Molecular weight of each segment: 
Segment_mw = {
    "PMMA": 100.1,
    "PS": 104.1,
    "PB": 54.1
}

# C_{i} per unit M_{w} 
c_i_coeff = {
    "PMMA": 0.0092,
    "PS": 0.0066,
}

# Molar group volumes of all relevant groups
R_k = {
    "CH=CH":1.1167,
    "CH2":0.6744,
    "C":0.2195,
    "CH3COO":1.9031,
    "CH3":0.9011,
    "ACH":0.5313,
    "ACCH":0.8121
}

# Number of each group in the polymer repeating unit
Polymer_groups = {
    "PB":{
        "CH=CH":1,
        "CH2":2
    },
    "PS":{
        "CH3":1,
        "CH2":1,
        "C":1,
        "CH3COO":1
    },
    "PMMA":{
        "ACH": 5,
        "ACCH":1,
        "CH2":1
    }
}

# Setting up volume expansion dicts
# rho_reference has units of g/ cm^3
rho_reference = {
    "PMMA":1.190,
    "PS":1.05,
    "PB":0.915
}

#units of Celcius
rho_reference_temp = {
    "PMMA":20.0,
    "PS":20.0,
    "PB":25.0
}

# Glass transition temperatures
T_G = {
    "PMMA": 105,
    "PS": 80
}

# volume thermal expansion coefficients: 
expansion_coefficients = {
    "PMMA":{
        "LEQ_TG":2.60e-4,
        "GEQ_TG":5.60e-4
    },
    "PS":{
        "LEQ_TG":1.7e-4,
        "GEQ_TG":5.1e-4
    },
    "PB":6.99e-4
}

"""
This section deals with calculating the combinatorial contribution

This requires calculating \phi_{i} which itself is dependent on r_{i}
"""
# The input to this function is the str of the species of interest. 
# It will then use this as the key to obtain the dictionary values from the necessary dicts
def r_i (Species):
    r = 0
    for group in list(Polymer_groups[Species].keys()):
        r = r + Polymer_groups[Species][group]*R_k[group]*Polymer_length[Species]
    
    return r

print (r_i("PMMA"))

def phi_i (x_i, r_i_1, r_i_2):
    phi = (x_i*(r_i_1**(0.75)))/((x_i*(r_i_1**(0.75)))+((1.0-x_i)*(r_i_2**(0.75))))

    return phi

def ln_gamma_combi_i (phi, x_i):
    ln_gamma_combi = ln(phi / x_i) + 1 - (phi / x_i)

    return ln_gamma_combi
    

"""
This section deals with the free volume contribution i.e. \gamma_{FV}

We need to calculate C_{i}, \tilde{V}_{i} and \tilde{V}_{m}
"""

# V_{i} has units of cm^/g 
def volume_i (Species, Temp):
    if Species in list(T_G.keys()):
        if Temp <= T_G[Species]:
            coeff = expansion_coefficients[Species]["LEQ_TG"]
        else:
            coeff = expansion_coefficients[Species]["GEQ_TG"]
    
    else:
        coeff = expansion_coefficients[Species]
    
    rho = rho_reference[Species] / (1 + coeff*(Temp - rho_reference_temp[Species]))

    volume = 1 / rho

    return volume

# Could possibly refactor red_volume_i to generalise the inputs to simply those as defined in red_volume_mix
def red_volume_i (vol, r_i_1):
    red_volume = vol/(15.17*1.28*r_i_1)

    return red_volume

def red_volume_mix (species_1, species_2, Temp, x_1):
    # Calculate weight fractions from mole fractions
    w_1 = x_1*(Polymer_length[species_1]*Segment_mw[species_1])/((x_1*Polymer_length[species_1]*Segment_mw[species_1])+((1.0-x_1)*Polymer_length[species_2]*Segment_mw[species_2]))
    w_2 = (1.0-x_1)*(Polymer_length[species_2]*Segment_mw[species_2])/((x_1*Polymer_length[species_1]*Segment_mw[species_1])+((1.0-x_1)*Polymer_length[species_2]*Segment_mw[species_2]))
    
    red_vol_mix = (volume_i(species_1, Temp)*w_1 + volume_i(species_2, Temp)*w_2)/(15.17*1.28*(r_i(species_1)*w_1 + r_i(species_2)*w_2))

    return red_vol_mix

def C_i_coeff (species):
    if species == "PMMA":
        C_1 = c_i_coeff[species]*Polymer_length[species]*Segment_mw[species]
    
    elif species == "PS":
        C_1 = c_i_coeff[species]*Polymer_length[species]*Segment_mw[species]
    
    elif species == "PB":
        C_1 = -0.640 + ((0.6744*0.146*2.0) + (1.1167*0.304))*Polymer_length[species]

    return C_1

def ln_gamma_fv_i(species_1, species_2, Temp, x_1):

    red_v_i_third = (red_volume_i(volume_i(species_1, Temp), r_i(species_1)))**(1.0/3.0)

    red_v_mix_third = (red_volume_mix(species_1, species_2, Temp, x_1))**(1.0/3.0)

    red_v_i_ = (red_volume_i(volume_i(species_1, Temp), r_i(species_1)))

    red_v_mix_ = (red_volume_mix(species_1, species_2, Temp, x_1))

    ln_gamma_fv = 3.0*C_i_coeff(species_1) * ln((red_v_i_third- 1.0)/(red_v_mix_third - 1.0)) - C_i_coeff(species_1)*((red_v_i_/red_v_mix_)-1.0)*((1.0 - (1.0/red_v_i_third))**(-1.0))
    
    return ln_gamma_fv

print(ln_gamma_fv_i("PS", "PB", 30.0, 0.5))


"""
This section is to evaluate the residual contribution to the activity coefficient
"""





"""
This section is to add up the respective contributions from the activity coeffcient 

"""

