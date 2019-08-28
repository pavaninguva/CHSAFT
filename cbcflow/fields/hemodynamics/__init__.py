"""
from AWSS import AWSS
from ICI import ICI
from LNWSS import LNWSS
from LSA import LSA
from MWSS import MWSS
from AOSI import AOSI
from PLC import PLC
from SCI import SCI
from VDR import VDR
"""
hemodynamic_fields = ["AWSS", "ICI", "LNWSS", "LSA", "MWSS", "AOSI", "PLC", "SCI", "VDR", "MinWSS"]

for f in hemodynamic_fields:
    exec("from cbcflow.fields.hemodynamics.%s import %s" %(f, f))
