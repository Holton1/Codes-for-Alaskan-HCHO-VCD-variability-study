import requests
import pandas as pd
import numpy as np
import datetime
import glob
import matplotlib.pyplot as plt
import matplotlib as mpl
import json
import os
import xarray as xr
import cmaps
import cmocean
import xesmf as xe
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

import regionmask
# mask out land region
land = regionmask.defined_regions.natural_earth.land_110

import netCDF4 as nc4
from netCDF4 import Dataset

import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

mpl.rcParams['font.sans-serif'] = "Arial"
mpl.rcParams['font.family'] = "sans-serif"
plt.rcParams['figure.figsize'] = (10, 6)
mpl.rc('xtick', labelsize=18) 
mpl.rc('ytick', labelsize=18)
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.titlesize'] = 18

figsize=(10,6)


#################################################################################




#############################################
#.  load ATom VOCs
#############################################

def load_ATom_VOCs_AK():
    
    #---------------------------------------------------
    #                   load HCHO
    #---------------------------------------------------
    fdir_ATom = "/import/GREENING/tzhao/ATom/"
    # import HCHO and VOC measurements from ATom-1


    # read in HCHO data (ISAF)
    fname = "MER-1HZ_DC8_ATom-1.nc"

    fh = Dataset(fdir_ATom + fname , mode='r')
    ATom_HCHO = np.array( fh.groups['ISAF-H2CO'].variables['CH2O'][:] )
    ATom_time_HCHO = np.array( fh.variables['time'][:] ) + 1451606400
    ATom_lat_HCHO = np.array( fh.groups['MMS'].variables['G_LAT'][:] )
    ATom_lon_HCHO = np.array( fh.groups['MMS'].variables['G_LONG'][:] )
    ATom_alt_HCHO = np.array( fh.groups['MMS'].variables['G_ALT'][:] )   # unit: m,   _FillValue: -99999
    ATom_pres_HCHO = np.array( fh.groups['MMS'].variables['P'][:] )   # unit: hPa,   _FillValue: -99999 
    ATom_temp_HCHO = np.array( fh.groups['MMS'].variables['T'][:] )


    ATom_HCHO[ATom_HCHO==-99999]=np.nan
    ATom_lat_HCHO[ATom_lat_HCHO==-99999]=np.nan
    ATom_lon_HCHO[ATom_lon_HCHO==-99999]=np.nan
    ATom_alt_HCHO[ATom_alt_HCHO==-99999]=np.nan
    ATom_pres_HCHO[ATom_pres_HCHO==-99999]=np.nan
    ATom_temp_HCHO[ATom_temp_HCHO==-99999]=np.nan

    # select region to AK
    x = np.where((ATom_lat_HCHO>=50)&(ATom_lat_HCHO<=75))
    y = np.where((ATom_lon_HCHO>=-170)&(ATom_lon_HCHO<=-130))
    iAK = np.intersect1d(x,y)

    ATom_HCHO = ATom_HCHO[iAK]
    ATom_time_HCHO = ATom_time_HCHO[iAK]
    ATom_lat_HCHO = ATom_lat_HCHO[iAK]
    ATom_lon_HCHO = ATom_lon_HCHO[iAK]
    ATom_alt_HCHO = ATom_alt_HCHO[iAK]
    ATom_pres_HCHO = ATom_pres_HCHO[iAK]
    ATom_temp_HCHO = ATom_temp_HCHO[iAK]
    
    
    # Pack in DataFrame

    # datetime obj
    datetime_obj = list(map(lambda a: datetime.datetime.utcfromtimestamp(a),  ATom_time_HCHO))

    # merge variables in df
    df = pd.concat(  [
        pd.DataFrame(datetime_obj, columns =['time']) ,\
        pd.DataFrame(np.array([ATom_lat_HCHO,ATom_lon_HCHO,ATom_alt_HCHO,ATom_pres_HCHO,ATom_temp_HCHO,ATom_HCHO]).T,\
                     columns =['lat','lon','alt','press','temp','HCHO_pptv'])
        ], axis=1, sort=False)
    ATom_HCHO_df = df.set_index(["time"])
    ATom_HCHO_df = ATom_HCHO_df.resample("1Min").mean()
    
    #---------------------------------------------------
    #                   load ISOP
    #---------------------------------------------------

    # read in ISOP data (WAS)
    fname = "MER-WAS_DC8_ATom-1.nc"

    fh = Dataset(fdir_ATom + fname , mode='r')
    ATom_ISOP = np.array( fh.groups['WAS'].variables['Isoprene_WAS'][:] )
    ATom_time_ISOP = np.array( fh.variables['time'][:] ) + 1451606400
    ATom_lat_ISOP = np.array( fh.groups['MMS'].variables['G_LAT'][:] )
    ATom_lon_ISOP = np.array( fh.groups['MMS'].variables['G_LONG'][:] )
    ATom_alt_ISOP = np.array( fh.groups['MMS'].variables['G_ALT'][:] )
    ATom_pres_ISOP = np.array( fh.groups['MMS'].variables['P'][:] )


    ATom_ISOP[ATom_ISOP==-99999]=np.nan
    ATom_lat_ISOP[ATom_lat_ISOP==-99999]=np.nan
    ATom_lon_ISOP[ATom_lon_ISOP==-99999]=np.nan
    ATom_alt_ISOP[ATom_alt_ISOP==-99999]=np.nan
    ATom_pres_ISOP[ATom_pres_ISOP==-99999]=np.nan

    ATom_ISOP[ATom_ISOP==-777]=np.nan
    ATom_lat_ISOP[ATom_lat_ISOP==-777]=np.nan
    ATom_lon_ISOP[ATom_lon_ISOP==-777]=np.nan
    ATom_alt_ISOP[ATom_alt_ISOP==-777]=np.nan
    ATom_pres_ISOP[ATom_pres_ISOP==-777]=np.nan

    ATom_ISOP[ATom_ISOP==-8888]=0.0


    # select region to AK
    x = np.where((ATom_lat_ISOP>=50)&(ATom_lat_ISOP<=75))
    y = np.where((ATom_lon_ISOP>=-170)&(ATom_lon_ISOP<=-130))
    iAK = np.intersect1d(x,y)

    ATom_ISOP = ATom_ISOP[iAK]
    ATom_time_ISOP = ATom_time_ISOP[iAK]
    ATom_lat_ISOP = ATom_lat_ISOP[iAK]
    ATom_lon_ISOP = ATom_lon_ISOP[iAK]
    ATom_alt_ISOP = ATom_alt_ISOP[iAK]
    ATom_pres_ISOP = ATom_pres_ISOP[iAK]
    
    
    # Pack in DataFrame

    # datetime obj
    datetime_obj = list(map(lambda a: datetime.datetime.utcfromtimestamp(a),  ATom_time_ISOP))

    # merge variables in df
    df = pd.concat(  [
        pd.DataFrame(datetime_obj, columns =['time']) ,\
        pd.DataFrame(np.array([ATom_lat_ISOP,ATom_lon_ISOP,ATom_pres_ISOP,ATom_ISOP]).T,\
                     columns =['lat','lon','press','ISOP_pptv'])
        ], axis=1, sort=False)
    ATom_ISOP_df = df.set_index(["time"])
    ATom_ISOP_df = ATom_ISOP_df.resample("1Min").mean()

    

    #---------------------------------------------------
    #                   load Monoterpenes
    #---------------------------------------------------        
    
    # read in Monoterpenes data (WAS)

    fname = "pinenes.csv"
    fh = pd.read_csv(fdir_ATom + fname)
    ATom_MTP = fh[["alpha-pinene","beta-pinene"]]
    ATom_time_MTP = (fh["UTC_Start_WAS"] + fh["UTC_Stop_WAS"])/2 + 1470009600

    ATom_MTP[ATom_MTP==-888]=0.0
    ATom_MTP[ATom_MTP==-99999]=np.nan
    ATom_MTP[ATom_MTP==-777]=np.nan
    
    
    # Pack in DataFrame

    # datetime obj
    datetime_obj = list(map(lambda a: datetime.datetime.utcfromtimestamp(a),  ATom_time_MTP))

    # merge variables in df
    df = pd.concat(  [
        pd.DataFrame(datetime_obj, columns =['time']) ,\
        ATom_MTP,\
        ], axis=1, sort=False)
    ATom_MTP_df = df.set_index(["time"])
    ATom_MTP_df = ATom_MTP_df.resample("1Min").mean()
    
    
    #------------------------------------------------------------
    #            load MVK and MACR from TOGA measurements
    #------------------------------------------------------------      
    
    # read in MVK and MACR data (TOGA)
    filelist = sorted( glob.glob(fdir_ATom + "TOGA/" + "*"))
    frames = []
    frames_time = []

    for ifile in filelist:
        fh = pd.read_csv( ifile , skiprows=94)
        ATom_VOCs = fh[["MVK_TOGA","MAC_TOGA","CH2O_TOGA","Isoprene_TOGA","aPinene_TOGA","bPineneMyrcene_TOGA"]]
        ATom_time_VOCs = (fh["Time_Start"] + fh["Time_Stop"])/2 + 1470009600

        ATom_VOCs[ATom_VOCs==-999]=np.nan
        ATom_VOCs[ATom_VOCs==-777]=np.nan
        ATom_VOCs[ATom_VOCs==-888]=0.0
    
        frames.append(ATom_VOCs)
        frames_time.append(ATom_time_VOCs)

    ATom_VOCs_all = pd.concat(frames, axis=0).reset_index(drop=True)
    ATom_time_VOCs_all = pd.concat(frames_time, axis=0).reset_index(drop=True)
    
    
    # Pack in DataFrame

    # datetime obj
    datetime_obj = list(map(lambda a: datetime.datetime.utcfromtimestamp(a),  ATom_time_VOCs_all))

    # merge variables in df
    df = pd.concat(  [
        pd.DataFrame(datetime_obj, columns =['time']) ,\
        ATom_VOCs_all,\
        ], axis=1, sort=False)
    ATom_VOCs_all_df = df.set_index(["time"])
    ATom_VOCs_all_df = ATom_VOCs_all_df.resample("1Min").mean()
    
    
    
    
    
    #---------------------------------------------------
    #             merge all VOCs in one Df
    #---------------------------------------------------  

    frames = [ATom_HCHO_df,ATom_ISOP_df["ISOP_pptv"],ATom_MTP_df, ATom_VOCs_all_df]
    ATom_VOC_df = pd.concat(frames, axis=1, sort=False)
    ATom_VOC_df[ ~np.isnan(ATom_VOC_df["lat"]) ]
    ATom_VOC_df = ATom_VOC_df.reset_index()
    
    # convert to xarray
    ATom_VOC_da = ATom_VOC_df.to_xarray()
    tmp = ATom_VOC_da.assign_coords(time=ATom_VOC_da.time)
    tmp = tmp.sortby(tmp.time)
    ATom_VOC_da = tmp.drop_vars('index')
    
    return ATom_VOC_da


# ##############   read data   ###############


# # load from rum
# ATom_VOC_da = rum.ATom.load_ATom_VOCs_AK()



    


    
def load_GC_ATom_VOCs_AK():    

    fdir_ATom_Planeflight = "/import/GREENING/tzhao/planeflight_2016/nest08_merra2_2x25_tropchem/"
    # load GEOS-Chem planelog

    filelist = sorted( glob.glob( fdir_ATom_Planeflight + "plane.log.2016080*" ))
    frames = []
    frames_time = np.array([])

    for ifile in filelist:
        fh = pd.read_csv( ifile ,sep='\s+' )
        ATom_GC = fh[["CH2O","ISOP","MTPA","MVK","MACR","GMAO_TEMP","GMAO_PRES","LAT","LON","P-I"]]

        YYYYMMDDHHMM = [i+j for i, j in zip( map(str, fh["YYYYMMDD"]), map(lambda x:str(x).zfill(4), fh["HHMM"]) )]
        ATom_time_GC = np.array( list( map( lambda x: datetime.datetime.strptime(x, '%Y%m%d%H%M'), YYYYMMDDHHMM ) )  )

    #     ATom_GC[ATom_GC<=-777]=np.nan
        frames.append(ATom_GC)
        frames_time = np.append(frames_time, ATom_time_GC)

    ATom_GC_all = pd.concat(frames, axis=0).reset_index(drop=True)





    # Convert GEOS-Chem planelog to netcdf

    # merge variables in df
    df = pd.concat(  [
        pd.DataFrame(frames_time, columns =['time']) ,\
        ATom_GC_all,\
        ], axis=1, sort=False)
    ATom_GC_all_df = df.set_index(["time"])
    ATom_GC_all_df = ATom_GC_all_df.resample("1Min").mean()


    # convert to xarray
    ATom_GC_all_da = ATom_GC_all_df.to_xarray()
    tmp = ATom_GC_all_da.assign_coords(time=ATom_GC_all_da.time)
    ATom_GC_all_da = tmp.sortby(tmp.time)
    
    return ATom_GC_all_da