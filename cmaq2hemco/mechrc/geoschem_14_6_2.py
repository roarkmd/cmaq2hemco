__all__ = ['writeconfig']


cq2gc = {
    'ACET': [['ACET', '']],
    'ACR': [['ACR', '']],
    'ACTA': [['ACTA', '']],
    'ALD2': [['ALD2', '']],
    'ALK4': [['ALK4', '']],
    'ALK6': [['ALK6', '']],
    'BALD': [['BALD', '']],
    'BC': [['BCPI', '70'], ['BCPO', '71']],
    'BENZ': [['BENZ', '']],
    'C2H2': [['C2H2', '']],
    'C2H4': [['C2H4', '']],
    'C2H6': [['C2H6', '']],
    'C3H8': [['C3H8', '']],
    'C4H6': [['C4H6', '']],
    'CCl4': [['CCl4', '']],
    'CH2Br2': [['CH2Br2', '']],
    'CH2Cl2': [['CH2Cl2', '']],
    'CH2O': [['CH2O', '']],
    'CH3Br': [['CH3Br', '']],
    'CH3Cl': [['CH3Cl', '']],
    'CH3I': [['CH3I', '']],
    'CH4': [['CH4', '']],
    'CHCl3': [['CHCl3', '']],
    'CL2': [['CL2', '']],
    'CO': [['CO', '']],
    'CSL': [['CSL', '']],
    'EBZ': [['EBZ', '']],
    'EOH': [['EOH', '']],
    'FURA': [['FURA', '']],
    'GLYX': [['GLYX', '']],
    'HACTA': [['HACTA', '']],
    'HCL': [['HCL', '']],
    'HCOOH': [['HCOOH', '']],
    'HONO': [['HONO', '']],
    'ISOP': [['ISOP', '']],
    'IVOC': [['IVOC', '']],
    'LIMO': [['LIMO', '']],
    'MACR': [['MACR', '']],
    'MEK': [['MEK', '']],
    'MGLY': [['MGLY', '']],
    'MOH': [['MOH', '']],
    'MTPA': [['MTPA', '']],
    'MTPO': [['MTPO', '']],
    'MVK': [['MVK', '']],
    'NAP': [['NAP', '']],
    'NH3': [['NH3', '']],
    'NH4': [['NH4', '']],
    'NIT': [['NIT', '']],
    'NO': [['NO', '115']],
    'NO2': [['NO2', '']],
    'OC': [['OCPI', '72'], ['OCPO', '73']],
    'OCS': [['OCS', '']],
    'pFe': [['pFe', '']],
    'PHEN': [['PHEN', '']],
    'PNA': [['PNA', '']],
    'PRPE': [['PRPE', '']],
    'RCHO': [['RCHO', '']],
    'RCOOH': [['RCOOH', '']],
    'ROH': [['ROH', '']],
    'SO2': [['SO2', '']],
    'SO4': [['SO4', '']],
    'STYR': [['STYR', '']],
    'SULF': [['SULF', '']],
    'TMB': [['TMB', '']],
    'TOLU': [['TOLU', '']],
    'XYLE': [['XYLE', '']],
}
# ignore special species, inventory meta variables, HAP tracers, UNK/UNR
for key in ['NH3_FERT','HFLUX', 'VOC_INV','NMOG','CH4_INV','ACROLEIN','BUTADIENE13','ETHYLBENZ','UNK','UNR']:
    cq2gc[key] = []

# Ignore unused particulate species
for key in ['PAL', 'PCA', 'PCL', 'PH2O', 'PK', 'PSI', 'PTI','PMC', 'PMG', 'PMN', 'PMOTHR', 'PNCOM']:
    cq2gc[key] = []

def writeconfig(outpath, year, sector, filepatt):
    from .core import writeconfig as genwriteconfig
    return genwriteconfig(outpath, year, sector, filepatt, cq2gc)
