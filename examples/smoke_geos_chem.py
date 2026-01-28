import pandas as pd
import numpy as np
import os
from cmaq2hemco.utils import pt2hemco, gd2hemco, gd2hemco_fast, gd_file, open_file
from cmaq2hemco.mechrc.geoschem_14_6_2 import writeconfig

def geos_chem_mw():
    molecDct = {"CO": 28, "NH3": 17, "NO": 46, "NO2": 46, "SO2": 64, \
          "SULF": 98, "BC": 1, "OC": 1, "PMFINE": 1, "NIT": 1, \
          "SO4": 1, "PMC": 1, "HGIIGAS": 200.59,\
          "HGNRVA":200.59, "PHGI": 1, "DIESEL_PMEC": 1, \
          "DIESEL_PMFINE": 1, "DIESEL_PMNO3":1, "DIESEL_PMOC": 1, \
          "DIESEL_PMSO4": 1, "DIESEL_PMC": 1, "ACROLEIN": 56.0633, \
          "BUTADIENE13":54.0904, "VOC_INV": 1.0, \
          "CL2": 70.91,"CHROMHEX_C": 1, \
          "CHROMHEX_F": 1, "CHROMTRI_C": 1, "CHROMTRI_F": 1,\
          "NICKEL_C": 1, "NICKEL_F": 1, "BERYLLIUM_C": 1, \
          "BERYLLIUM_F": 1,"CADMIUM_C": 1, "CADMIUM_F": 1, "LEAD_C":1,\
          "LEAD_F": 1, "MANGANESE_C":1, "MANGANESE_F": 1, \
          "OXYL": 106.165, "PXYL": 106.165, "MXYL": 106.165,\
           "TOLU": 92.1384, "CL4_ETHE": 165.83, "TRIETHYLAMINE": 101.19,\
           "HEXAMETHY_DIIS": 168.2, "CHCL3": 119.3776, "CL_ETHE": 62.5,\
           "CL4_ETHANE1122": 167.85, "ETOX": 44.0526, "QUINOLINE": 129.16,\
           "ACRYLONITRILE": 53.06, "CL2_C2_12": 98.9592, \
           "BR2_C2_12": 187.86,"HYDRAZINE": 32.05, "CARBONTET": 153.8227,\
           "DICHLOROPROPENE": 110.97,"PROPDICHLORIDE": 112.9, \
           "MAL_ANHYDRIDE": 98.06, "DICHLOROBENZENE":147.0002, \
           "TOL_DIIS": 174.1561, "CL2_ME": 84.93, "CL3_ETHE": 131.3883,\
           "HCL": 36.46, "HONO": 46, "NOX": 46, "PM2_5": 1, "PM10": 1, "HFLUX": 1,\
           "NH3_FERT": 17, "PAL": 1, "PCA": 1, "PCL": 1, \
           "pFe": 1, "PH2O": 1, "PK": 1, "PMG": 1, "PMN": 1, "PMOTHR": 1, \
           "PNA": 1, "PNCOM": 1, "NH4": 1, "PSI": 1, "PTI": 1,
           "ARSENIC_C": 1, "ARSENIC_F": 1, \
           "PAH_000E0": 379.00, "PAH_101E2": 268.00, "PAH_114E1": 256.00,
           "PAH_176E2": 302.00, "PAH_176E3": 244.00, "PAH_176E4": 248.00,
           "PAH_176E5": 228.00, "PAH_192E3": 278.00, "PAH_880E5": 196.00,
           "ACETONITRILE": 41.05, "ACRYLICACID": 72.06, "ACRYLONITRILE": 53.06,
           "CARBSULFIDE": 60.07, "CHLOROPRENE": 88.54, "ETHYLBENZ": 106.165,
           "HEXANE": 86.175, "METHCHLORIDE": 50.49, "STYRENE": 104.15,
           "XYLENES": 106.165, "NOX_INV": 46, "HF": 20.01, "NMOG": 1.0,
           "VOC_BEIS": 46, "APIN": 136.234, "BPIN": 120, "SESQ": 180, "NR": 24,
           "BENZOAPYRNE": 252.316, "TOG_INV": 1, "CO2": 44.01, "CO2_INV": 1,
           "N2O_INV": 1, "CH4_INV": 16.042,
           "ACET": 58.09, "ACR": 56.06, "ACTA": 60.06, "ALD2": 44.06, "ALK4": 58.12,
           "ALK6": 100.2, "BALD": 106.12, "BENZ": 78.12, "C2H2": 26.05, "C2H4": 28.05,
           "C2H6": 30.08, "C3H8": 44.11, "C4H6": 54.09, "CCl4": 153.82, "CH2Br2": 173.83,
           "CH2Cl2": 84.93, "CH2O": 30.03, "CH4": 16.04, "CHCl3": 119.35, "CSL": 108.14,
           "EBZ": 106.167, "EOH": 46.07, "FURA": 68.07, "GLYX": 58.04, "HCOOH": 46.03,
           "ISOP": 68.13, "IVOC": 240.5, "MEK": 72.11, "MOH": 32.05, "MTPA": 136.26,
           "MTPO": 136.26, "MVK": 70.09, "NAP": 128.18, "OCS": 60.07, "PHEN": 94.11,
           "PRPE": 42.09, "RCHO": 58.09, "RCOOH": 74.09, "ROH": 60.11, "STYR": 104.1491,
           "TMB": 106.167, "TOLU": 92.15, "UNK": 142.3, "UNR": 142.3, "XYLE": 106.18,
           "CH3Br": 94.94, "CH3Cl": 50.45, "CH3I": 141.94, "HACTA": 76.0514, "LIMO": 136.26,
           "MACR": 70.1, "MGLY": 72.07}
    df = pd.DataFrame({"SPECIES": molecDct.keys(), "MOLWT": molecDct.values()})
    df.to_csv("cmaq_geoschem14_6_3_molwt.csv", index=False)



def main():
    # Debug set to True will only process a single date
    # or whatever is in the debug if satement below
    debug = True

    # excluding burn and biogenics, and point fire sectors
    # openburn
    # beis4
    gkeys = """rwc
    rail
    onroad_gas
    onroad_diesel
    onroad_ca_adj_gas
    onroad_ca_adj_diesel
    np_solvents
    np_oilgas
    nonroad_gas
    nonroad_diesel
    nonpt
    mexico_onroad
    livestock
    canada_ptdust_adj
    canada_onroad
    canada_og2D
    canada_afdust_adj
    airports
    afdust_adj""".split()

    pkeys = """ptegu
    canmex_point
    cmv_c1c2_12
    cmv_c3_12
    pt_oilgas
    ptnonipm""".split()

    #Generate speciation molecular weights
    geos_chem_mw()

    # Date range to process emissions
    dates = pd.date_range('2022-01-01', '2022-12-31', freq='d')

    # Define grid by edges
    # Grid is latitude/longitude grid tenth degree grid
    elat = np.linspace(15, 65, 501)
    elon = np.linspace(-135, -50, 851)

    if debug:
        dates = pd.to_datetime(['2022-01-04'])
        # pkeys = pkeys[:1]
        # gkeys = gkeys[:1]
        print('**WARNING: in debug mode; only processing')
        print(dates)
        print(pkeys)
        print(gkeys)

    for date in dates:
        for gkey in gkeys:
            outpath = (
                f'epa2022v2/{gkey}/{gkey}_{date:%Y-%m-%d}_geoschem14_6_3_2022he_geoschem.ncf'
            )
            if os.path.exists(outpath):
                continue
            print(date, gkey)
            os.makedirs(os.path.dirname(outpath), exist_ok=True)
            try:
                # open_gdemis could fail if date is not available (eg., only weekday)
                gf = gd_file(open_file(date, tmpl='/home/mroark/git_repositories/cmaq2hemco/examples/rail/emis_mole_rail_%Y%m%d_12US1_geoschem14_6_3_2022he_geoschem.ncf'))
            except Exception as e:
                print(f'**WARNING:: Skipping {date} {gkey}: {e}')
                continue
            # using bilinear interpolation of fluxes
            rgf = gd2hemco_fast(outpath, gf, elat, elon, gc='geoschem14_6_3', nr='geoschem14_6_3')
            # use matrix interpolation for fractional area overlap (slow)
            # rgf = gd2hemco(outpath, gf, elat, elon, matrix=matrix)
            del rgf, gf
        
        for pkey in pkeys:
            outpath = (
                f'epa2022v2/{pkey}/{pkey}_{date:%Y-%m-%d}_geoschem14_6_3_2022he_geoschem.ncf'
            )
            if os.path.exists(outpath):
                continue
            print(date, pkey)
            os.makedirs(os.path.dirname(outpath), exist_ok=True)
            try:
                # open_ptemis could fail if date is not available (eg., only weekday)
                pf = gd_file(open_file(date, tmpl='/home/mroark/git_repositories/cmaq2hemco/examples/rail/emis_mole_rail_%Y%m%d_12US1_geoschem14_6_3_2022he_geoschem.ncf'))
                #pf = gd_file(open_date(date, tmpl='sector/emis_mole_rail_%Y%m%d_12US1_geoschem14_6_3_2022he_geoschem.ncf', \
                #                       bucket='/home/mroark/git_repositories/cmaq2hemco/examples'))
            except IOError as e:
                print(f'**WARNING:: Skipping {date} {pkey}: {e}')
                continue
            rpf = pt2hemco(outpath, pf, elat, elon, gc='geoschem14_6_3', nr='geoschem14_6_3')  # apply plume rise
            del rpf, pf

    for sector in gkeys + pkeys:
        hcpath = f'epa2022v2/{sector}/HEMCO_{sector}.rc'
        secttmpl = f'epa2022v2/{sector}/{sector}_%Y-%m-%d_geoschem14_6_3_2022he_geoschem.ncf'
        writeconfig(hcpath, 2022, sector, secttmpl)



if __name__ == '__main__':
    main()
