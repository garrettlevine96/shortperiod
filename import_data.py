import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import ascii

'''
SET OF FUNCTIONS TO IMPORT EXOPLANET DATA FROM
A FEW DIFFERENT SOURCES FOR USE IN THE PEAS-
IN-A-POD PROJECT.
'''

def rename_planets_to_Kepler(pl):
    pl[np.where(pl == 'K00841.03')] = 'Kepler-27'
    pl[np.where(pl == 'K00115.03')] = 'Kepler-105'
    pl[np.where(pl == 'K01563.03')] = 'Kepler-305'
    pl[np.where(pl == 'K01831.03')] = 'Kepler-324'
    pl[np.where(pl == 'K00332.02')] = 'Kepler-526'
    pl[np.where(pl == 'K00427.01')] = 'Kepler-549'
    pl[np.where(pl == 'K02162.02')] = 'Kepler-1126'
    return pl

def import_hadden(filename = 'data/hadden_2017.txt'):
    table = ascii.read(filename)
    planets = np.array(table['Planet'])
    planets = rename_planets_to_Kepler(planets)

    # drop the planet "b", "c", "d" for just the host name
    system = []
    for pl in planets:
        host = pl.split(" ")[0]
        system.append(host)

    # import the "robust" flag
    robust = np.array(table['robust'])

    # period [days]
    P = np.array(table['Per'])

    # planet radius [Earth radii]
    R_pl = np.array(table['Rad'])
    R_pl_err_low = np.array(table['e_Rad'])
    R_pl_err_high = np.array(table['E_Rad'])

    # planet mass, default prior [Earth masses]
    M_pl = np.array(table['Md'])
    M_pl_err_low = np.array(table['e_Md'])
    M_pl_err_high = np.array(table['E_Md'])

    # planet density, default prior [g/cm^3]
    rho_pl = np.array(table['rhod'])
    rho_pl_err_low = np.array(table['e_rhod'])
    rho_pl_err_high = np.array(table['E_rhod'])

    # planet mass, high prior [Earth masses]
    M_pl_high = np.array(table['Mh'])
    M_pl_high_err_low = np.array(table['e_Mh'])
    M_pl_high_err_high = np.array(table['E_Mh'])

    # planet density, high prior [Earth masses]
    rho_pl_high = np.array(table['rhoh'])
    rho_pl_high_err_low = np.array(table['e_rhoh'])
    rho_pl_high_err_high = np.array(table['E_rhoh'])

    # star mass [Solar masses]
    M_st = np.array(table['Mstar'])
    M_st_err_low = np.array(table['e_Mstar'])
    M_st_err_high = np.array(table['E_Mstar'])

    # initialize a pandas DataFrame
    cols = ["Name", "Host", "Robust", "P", "R_pl", "R_pl_err_low", "R_pl_err_high",\
        "M_pl", "M_pl_err_low", "M_pl_err_high", "rho_pl", "rho_pl_err_low",\
        "rho_pl_err_high", "M_pl_high", "M_pl_high_err_low", "M_pl_high_err_high",\
        "rho_pl_high", "rho_pl_high_err_low", "rho_pl_high_err_high",\
        "M_st", "M_st_err_low", "M_st_err_high"]
    df = pd.DataFrame(columns = cols)

    # add the rows to the DataFrame
    for i in range(len(planets)):
        listy = [planets[i], system[i], robust[i], P[i], R_pl[i], R_pl_err_low[i], R_pl_err_high[i],\
            M_pl[i], M_pl_err_low[i], M_pl_err_high[i], rho_pl[i], rho_pl_err_low[i],\
            rho_pl_err_high[i], M_pl_high[i], M_pl_high_err_low[i], M_pl_high_err_high[i],\
            rho_pl_high[i], rho_pl_high_err_low[i], rho_pl_high_err_high[i],\
            M_st[i], M_st_err_low[i], M_st_err_high[i]]
        df.loc[len(df)] = listy

    return df


def import_weiss(filename = "data/weiss_2018.txt"):
    """
    function to import the machine-readable data table
    from Weiss et al. (2018) "peas-in-a-pod" paper

    returns a pandas DataFrame object
    """
    # load the ASCII text file
    cks = np.loadtxt(filename, skiprows=30, dtype=object, delimiter = "\n")
    # define the column names
    cols = ["KOI_st", "KIC", "KOI_pl","Kepname", "Mass_st", "Radius_st",\
        "CDPP", "b", "Period_pl", "Radius_pl", "Radius_pl_sigma", "T_pl"]
    # initialize the pandas DataFrame
    weiss_df = pd.DataFrame(columns = cols)
    
    # populate the dataframe from the numpy array
    for i in range(len(cks)):
        # grab the row corresponding to the next planet
        row = cks[i]
        
        # grab the relevant fields
        # Kepler Object of Interest number of star
        koi = str(row[0:6]).strip()
        # Kepler Input Catalog number of star
        kic_st = str(row[7:15]).strip()
        # Kepler Object of Interest number of planet
        koi_pl = str(row[16:25]).strip()
        # Kepler name of planet, if named
        kepname = str(row[26:39]).strip()
        # Stellar mass (solar masses)
        mass_st = float(row[40:46])
        # Stellar radius (solar radii)
        radius_st = float(row[47:52])
        # Six-Hour Combined Differential 
        # Photometric Precision (ppm)
        cdpp = float(row[53:60])
        # Impact parameter
        b = float(row[61:67])
        # Orbital period
        period_pl = float(row[68:82])
        # Planet radius (Earth radii)
        radius_pl = float(row[83:98])
        # Uncertainty in planet radius
        radius_pl_sigma = float(row[99:114])
        # Equilibrium temperature [K]
        teq = float(row[115:129])
        
        # put these fields in a list
        listy = [koi, kic_st, koi_pl, kepname, mass_st, radius_st,\
            cdpp, b, period_pl, radius_pl, radius_pl_sigma, teq]
        
        # append this list as the next row in the pandas DataFrame
        weiss_df.loc[len(weiss_df)] = listy
        
    # return the DataFrame of values from Weiss et al. (2018)
    return weiss_df

def import_reinhold(filename = "data/reinhold_2013.txt"):
    table = np.genfromtxt(filename, dtype = "float",\
        missing_values = "---", skip_header = 67)
    cols = ["KIC", "Prot", "sigma_Prot",\
        "Prot2", "sigma_Prot2", "Teff",\
        "Logg", "Rvariability", "B-V"]
    df = pd.DataFrame(table, columns = cols)
    return df

def import_moons(filename = 'data/moons.csv'):
    moons = pd.read_csv(filename)
    moons = moons.drop([4, 5])
    return moons

def import_berger2020a(filename = "data/berger_2020a.txt"):
    cols = ["KIC", "KOI", "PD", "Radius", "E_Radius1", "E_Radius2",\
        "a", "E_a", "e_a", "Flux", "E_flux", "e_Flux", "ZAMSFlux", "Flag"]
    delims = (8, 8, 10, 9, 9, 9, 7, 7, 7, 10, 9, 9, 10, 14)
    df = pd.DataFrame(np.genfromtxt(filename,\
        skip_header = 32, dtype = "str", delimiter = delims))
    df.columns = cols 

    df['KIC'] = df["KIC"].astype(int)
    df['KOI'] = df['KOI'].astype(float)
    df['Radius'] = df["Radius"].astype(float)
    df['E_Radius1'] = df["E_Radius1"].astype(float)
    df['E_Radius2'] = df["E_Radius2"].astype(float)
    df['Flux'] = df["Flux"].astype(float)
    df['a'] = df['a'].astype(float)
    return df

def get_full_KOI_berger(array):
    full_KOIs = []
    for a in array:
        stringy = str(a).strip(" ")
        full_KOIs.append("K" + (8 - len(stringy)) * "0" + stringy)
    return np.array(full_KOIs)


def import_exoplanet_archive_KOIs(filename = "data/KOIs_archive_20211020.csv"):
    df = pd.read_csv(filename, skiprows = 47)
    df = df[["kepid", "kepoi_name", "kepler_name",\
        "koi_period", "koi_srad", "koi_srad_err1", "koi_srad_err2"]]
    
    return df


##############################


def fulton(P, m=-0.09, a=0.37):
    return 10**(m * np.log10(P) + a)

def super_Earth(P, R, m=-0.09, a=0.37):
    # find the maximum allowable radius
    # for a super-Earth
    maximum_SE = fulton(P, m=m, a=a)
    # compare the radius of the planet
    # to the maximum allowable radius
    return R < maximum_SE


