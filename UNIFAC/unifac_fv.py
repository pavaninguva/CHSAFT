"""
This section is for capturing the physical properties / group parameters of the polymers

We define the extent of polymerisation, number of groups etc
"""
# Update the input polymer properties as needed
Polymer_length = {
    "PMMA": 1000,
    "PS":1000,
    "PB":1000
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

def r_i (Species):
    r = 0
    for group in list(Polymer_groups[Species].keys()):
        r = r + Polymer_groups[Species][group]*R_k[group]*Polymer_length[Species]
    
    return r

print (r_i("PMMA"))

def phi_i ():
    var = ()



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

