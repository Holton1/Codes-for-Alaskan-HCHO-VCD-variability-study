import pandas as pd
import numpy as np
import xarray as xr
import os
import datetime
from datetime import datetime
from scipy import interpolate
from scipy.interpolate import griddata
import netCDF4 as nc4
from netCDF4 import Dataset
import xesmf as xe
import sys

def file_list(dirname, ext='.csv'):
    import os
    #获取目录下所有特定后缀的文件
    #@param dirname: str 目录的完整路径
    #@param ext: str 后缀名, 以点号开头
    #@return: list(str) 所有子文件名(不包含路径)组成的列表,
    return list(filter(
        lambda filename: os.path.splitext(filename)[1] == ext,
        os.listdir(dirname)))

import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.cm as cm
xlabelfont = {'family' : 'DejaVu Sans',
              'weight' : 'normal',
              'size'   : 16 ,
             }
vmin=0
vmax=1e16
cmap = mpl.cm.seismic
# norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

        
def calbox(clat,clon,km):
    import math

#     clat = 64.8
#     clon = -147.7A
#     km = 30

    cos_clat = math.cos(clat*math.pi/180)
    dlon = km/(111*cos_clat)
    dlat = km/111

    latN = clat + dlat/2
    latS = clat - dlat/2
    latrange = [latS, latN]

    lonE = clon + dlon/2
    lonW = clon - dlon/2
    lonrange = [lonW, lonE]
    
    return latrange, lonrange

selbox = lambda da: da.where( (da.lat>latrange[0]) & (da.lat<latrange[1]) &\
                             (da.lon>lonrange[0]) & (da.lon<lonrange[1]) ).mean(dim='lat').mean(dim='lon')
sellatlon = lambda da: da.where( (da.lat>latrange[0]) & (da.lat<latrange[1]) &\
                             (da.lon>lonrange[0]) & (da.lon<lonrange[1]) )



#--------------------------------------------------------------------------------
#                        Start the process
#--------------------------------------------------------------------------------

runname = sys.argv[1]
YYYY = sys.argv[2]


jndirSIF="/import/GREENING/tzhao/jndata/TROPOMI_SIF/"
jndir = "/import/GREENING/tzhao/jndata/TROPOMI_HCHO/"
jndirCO = "/import/GREENING/tzhao/jndata/TROPOMI_CO/"
jndirPdr = "/import/GREENING/tzhao/jndata/Pandora/"
jndirGC = "/import/GREENING/tzhao/jndata/GEOS-Chem/"
jndirGC_NoWF = jndirGC + "MEGANon_PFTol_NoWF/"
jndirGC_WF = jndirGC + "MEGANon_PFTol_WF/"

# outdir = "/import/GREENING/tzhao/jndata/output_Online_20200331/"
workdir = "/home/tzhao/my_jupyter_work/newest/"
# lat, lon 0.5x0.625
###############################################################################
lat_gc_AK = np.load(jndir+'lat_gc_AK_050625.npy')
lon_gc_AK = np.load(jndir+'lon_gc_AK_050625.npy')
# lon_tpm_AK = np.load(jndir+'lon_gc_AK_2x25Over16.npy')
# lat_tpm_AK = np.load(jndir+'lat_gc_AK_2x25Over16.npy')

TPM_hsigma = xr.open_dataset(jndir+"TROPOMI_hsigma.nc", engine="netcdf4")
b_tpm = TPM_hsigma["b_tpm"]
ap_tpm = TPM_hsigma["ap_tpm"]






if runname == "MEGANon_PFTol_WF":
    runname_dir = jndirGC_WF
elif runname == "MEGANon_PFTol_NoWF":
    runname_dir = jndirGC_NoWF


if YYYY=="2018":

    # Read GESO-Chem data
    # jndirGC_NoWF = jndirGC + "MEGANol_PFTol_NoWF/"

    SpeciesConc_CH2O = xr.open_dataarray(runname_dir+"SpeciesConc_CH2O_18.nc", engine="netcdf4").loc['2018-05-14T21:00:00':'2018-08-31T22:00:00']
    # SpeciesConc_ISOP = xr.open_dataarray(jndirGC_NoWF+"SpeciesConc_ISOP_18.nc", engine="netcdf4").loc['2018-05-14T21:00:00':'2018-08-31T22:00:00']

    Met_BXHEIGHT = xr.open_dataarray(jndirGC_WF+"Met_BXHEIGHT_18.nc", engine="netcdf4").loc['2018-05-14T21:00:00':'2018-08-31T22:00:00']
    Met_AIRDEN = xr.open_dataarray(jndirGC_WF+"Met_AIRDEN_18.nc", engine="netcdf4").loc['2018-05-14T21:00:00':'2018-08-31T22:00:00']
    Met_PBLH = xr.open_dataarray(jndirGC_WF+"Met_PBLH_18.nc", engine="netcdf4").loc['2018-05-14T21:00:00':'2018-08-31T22:00:00']
    Met_SLP = xr.open_dataarray(jndirGC_WF+"Met_SLP_18.nc", engine="netcdf4").loc['2018-05-14T21:00:00':'2018-08-31T22:00:00']
    Met_T = xr.open_dataarray(jndirGC_WF+"Met_T_18.nc", engine="netcdf4").loc['2018-05-14T21:00:00':'2018-08-31T22:00:00']
    ap_gc = xr.open_dataarray(jndirGC+"ap_gc.nc", engine="netcdf4")
    b_gc = xr.open_dataarray(jndirGC+"b_gc.nc", engine="netcdf4")


    # Read TROPOMI GCHiR dataset

    TPM_AvKapri_18 = xr.open_dataset(jndir+"TROPOMI_AvKapri_GCHiR_hr_18_Fullrange.nc", engine="netcdf4")
    TPM_18 = xr.open_dataset(jndir+"TROPOMI_HCHO_GCHiR_hr_18_Fullrange.nc", engine="netcdf4")


    AveragingKernel_TPM = TPM_AvKapri_18["AveragingKernel_TPM_GCHiR"].loc['2018-05-14T21:00:00':'2018-08-31T22:00:00']
    apriori_TPM = TPM_AvKapri_18["apriori_TPM_GCHiR"].loc['2018-05-14T21:00:00':'2018-08-31T22:00:00']

    HCHOVCD_TPM = TPM_18["HCHOVCD_TPM_GCHiR"].loc['2018-05-14T21:00:00':'2018-08-31T22:00:00']
    HCHOVCD_err_TPM = TPM_18["HCHOVCD_err_TPM_GCHiR"].loc['2018-05-14T21:00:00':'2018-08-31T22:00:00']
    cloudfrac_TPM = TPM_18["cloudfrac_TPM_GCHiR"].loc['2018-05-14T21:00:00':'2018-08-31T22:00:00']
    SZA_TPM = TPM_18["SZA_TPM_GCHiR"].loc['2018-05-14T21:00:00':'2018-08-31T22:00:00']
    VZA_TPM = TPM_18["VZA_TPM_GCHiR"].loc['2018-05-14T21:00:00':'2018-08-31T22:00:00']
    AMF_HCHO_trop_TPM = TPM_18["AMF_HCHO_trop_TPM_GCHiR"].loc['2018-05-14T21:00:00':'2018-08-31T22:00:00']
    qa_TPM = TPM_18["qa_TPM_GCHiR"].loc['2018-05-14T21:00:00':'2018-08-31T22:00:00']
    Psurf_TPM = TPM_18["Psurf_TPM_GCHiR"].loc['2018-05-14T21:00:00':'2018-08-31T22:00:00']/100
    albedo_TPM = TPM_18["albedo_TPM_GCHiR"].loc['2018-05-14T21:00:00':'2018-08-31T22:00:00']
    HCHOSCD_corr_TPM = TPM_18["HCHOSCD_corr_TPM_GCHiR"].loc['2018-05-14T21:00:00':'2018-08-31T22:00:00']
    HCHOVCD0_TPM = TPM_18["HCHOVCD0_TPM_GCHiR"].loc['2018-05-14T21:00:00':'2018-08-31T22:00:00']



    
elif YYYY == "2019":
    
    SpeciesConc_CH2O = xr.open_dataarray(runname_dir+"SpeciesConc_CH2O_19.nc", engine="netcdf4").loc['2019-05-01T00:00:00':'2019-08-31T23:00:00']
    # SpeciesConc_ISOP = xr.open_dataarray(jndirGC_NoWF+"SpeciesConc_ISOP_19.nc", engine="netcdf4").loc['2019-05-01T00:00:00':'2019-08-31T23:00:00']
    Met_BXHEIGHT = xr.open_dataarray(jndirGC_WF+"Met_BXHEIGHT_19.nc", engine="netcdf4").loc['2019-05-01T00:00:00':'2019-08-31T23:00:00']
    Met_AIRDEN = xr.open_dataarray(jndirGC_WF+"Met_AIRDEN_19.nc", engine="netcdf4").loc['2019-05-01T00:00:00':'2019-08-31T23:00:00']
    Met_PBLH = xr.open_dataarray(jndirGC_WF+"Met_PBLH_19.nc", engine="netcdf4").loc['2019-05-01T00:00:00':'2019-08-31T23:00:00']
    Met_SLP = xr.open_dataarray(jndirGC_WF+"Met_SLP_19.nc", engine="netcdf4").loc['2019-05-01T00:00:00':'2019-08-31T23:00:00']
    Met_T = xr.open_dataarray(jndirGC_WF+"Met_T_19.nc", engine="netcdf4").loc['2019-05-01T00:00:00':'2019-08-31T23:00:00']
    ap_gc = xr.open_dataarray(jndirGC+"ap_gc.nc", engine="netcdf4")
    b_gc = xr.open_dataarray(jndirGC+"b_gc.nc", engine="netcdf4")


    # Read TROPOMI GCHiR dataset

    TPM_AvKapri_19 = xr.open_dataset(jndir+"TROPOMI_AvKapri_GCHiR_hr_19_Fullrange__AddHiR.nc", engine="netcdf4")

    TPM_19 = xr.open_dataset(jndir+"TROPOMI_HCHO_GCHiR_hr_19_Fullrange__AddHiR.nc", engine="netcdf4")

    AveragingKernel_TPM = TPM_AvKapri_19["AveragingKernel_TPM_GCHiR"]
    apriori_TPM = TPM_AvKapri_19["apriori_TPM_GCHiR"]

    HCHOVCD_TPM = TPM_19["HCHOVCD_TPM_GCHiR"].loc['2019-05-01T00:00:00':'2019-08-31T23:00:00']
    HCHOVCD_err_TPM = TPM_19["HCHOVCD_err_TPM_GCHiR"].loc['2019-05-01T00:00:00':'2019-08-31T23:00:00']
    cloudfrac_TPM = TPM_19["cloudfrac_TPM_GCHiR"].loc['2019-05-01T00:00:00':'2019-08-31T23:00:00']
    SZA_TPM = TPM_19["SZA_TPM_GCHiR"].loc['2019-05-01T00:00:00':'2019-08-31T23:00:00']
    VZA_TPM = TPM_19["VZA_TPM_GCHiR"].loc['2019-05-01T00:00:00':'2019-08-31T23:00:00']
    AMF_HCHO_trop_TPM = TPM_19["AMF_HCHO_trop_TPM_GCHiR"].loc['2019-05-01T00:00:00':'2019-08-31T23:00:00']
    qa_TPM = TPM_19["qa_TPM_GCHiR"].loc['2019-05-01T00:00:00':'2019-08-31T23:00:00']
    Psurf_TPM = TPM_19["Psurf_TPM_GCHiR"].loc['2019-05-01T00:00:00':'2019-08-31T23:00:00']/100
    albedo_TPM = TPM_19["albedo_TPM_GCHiR"].loc['2019-05-01T00:00:00':'2019-08-31T23:00:00']
    HCHOSCD_corr_TPM = TPM_19["HCHOSCD_corr_TPM_GCHiR"].loc['2019-05-01T00:00:00':'2019-08-31T23:00:00']
    HCHOVCD0_TPM = TPM_19["HCHOVCD0_TPM_GCHiR"].loc['2019-05-01T00:00:00':'2019-08-31T23:00:00']
    
    

    

# GEOS-Chem noon time window

def GCnoon(var_da, hour1, hour2):
    var_da_noon = var_da.where( (var_da.time.dt.hour>=hour1)&(var_da.time.dt.hour<=hour2)  )
    return var_da_noon

# UTC window
hour1 = 12+8
hour2 = 15+8

# # GEOS-CHem HCHO VCD
# tmp_vcd = (SpeciesConc_CH2O*Met_AIRDEN*Met_BXHEIGHT).sum(dim='lev')*6.02e23*1e3/29/1e4
# HCHOVCD_GC = GCnoon(tmp_vcd, hour1, hour2)







# Recalculate TROPOMI AMF from GC profile
# TROPOMI
P_TPM = (Psurf_TPM*b_tpm + ap_tpm/100).transpose('time','lev','lat','lon')
# Pressure at each GC layers
P_center_GC = Met_SLP*b_gc + ap_gc/100
P_center_GC = np.array(P_center_GC).transpose((0,3,1,2))
P_center_GC = xr.DataArray(P_center_GC, coords=[Met_T.time, Met_T.lev, lat_gc_AK, lon_gc_AK], dims=['time','lev' ,'lat','lon'])





# Altitudes at each GC layers
H_edge_GC = np.zeros((Met_BXHEIGHT.time.size,Met_BXHEIGHT.lev.size+1,\
                      Met_BXHEIGHT.lat.size,Met_BXHEIGHT.lon.size))
H_center_GC = np.zeros(Met_BXHEIGHT.shape)
for ilev in range(0,Met_BXHEIGHT.lev.size):
    H_edge_GC[:,ilev+1,:,:] += Met_BXHEIGHT[:,ilev,:,:]
for ilev in range(0,Met_BXHEIGHT.lev.size):
    H_center_GC[:,ilev,:,:] = ( H_edge_GC[:,ilev,:,:] + H_edge_GC[:,ilev+1,:,:] )*0.5

H_center_GC = xr.DataArray(H_center_GC, coords=[Met_T.time, Met_T.lev, lat_gc_AK, lon_gc_AK], dims=['time','lev' ,'lat','lon'])






# T in TROPOMI levels
# BXH in TROPOMI levels
# GC Number Density in each TROPOMI levels

gc_ND_HCHO = SpeciesConc_CH2O*Met_AIRDEN*1e3/29*6.02e23
da_gc_ND_HCHO = xr.DataArray(gc_ND_HCHO, coords=[Met_T.time, Met_T.lev, lat_gc_AK, lon_gc_AK], dims=['time','lev' ,'lat','lon'])




T_TPM_empty = np.empty(P_TPM.shape, float)
BXHEIGHT_TPM_empty = np.empty(P_TPM.shape, float)
gc_ND_HCHO_TPMgrids_empty = np.empty(P_TPM.shape, float)


for it in range(0,Met_T.time.size):
    print("it = "+str(it)+"starts")
    if np.int32(P_TPM[it,:,:].sum())==0:
        continue
    for ilat in range(0,Met_T.lat.size):
        for ilon in range(0,Met_T.lon.size):
            
            P_TPM_onecolumn = P_TPM[it,:,ilat,ilon].drop_vars('time').drop_vars('lat').drop_vars('lon')
            if np.isnan(P_TPM_onecolumn).sum()>0:
                continue
            P_center_GC_onecolumn = P_center_GC[it,:,ilat,ilon].drop_vars('time').drop_vars('lat').drop_vars('lon')
            
            
            Met_T_onecolumn = Met_T[it,:,ilat,ilon].assign_coords(lev=P_center_GC_onecolumn)
            T_TPM_empty[it,:,ilat,ilon] = Met_T_onecolumn.interp(lev=P_TPM_onecolumn)
            Met_BXHEIGHT_onecolumn = Met_BXHEIGHT[it,:,ilat,ilon].assign_coords(lev=P_center_GC_onecolumn)
            BXHEIGHT_TPM_empty[it,:,ilat,ilon] = Met_BXHEIGHT_onecolumn.interp(lev=P_TPM_onecolumn)
            da_gc_ND_HCHO_onecolumn = da_gc_ND_HCHO[it,:,ilat,ilon].assign_coords(lev=P_center_GC_onecolumn)
            gc_ND_HCHO_TPMgrids_empty[it,:,ilat,ilon] = da_gc_ND_HCHO_onecolumn.interp(lev=P_TPM_onecolumn)


# remove 0
T_TPM_empty[T_TPM_empty==0]=np.nan
# pack in xarray
T_TPM = xr.DataArray(T_TPM_empty, coords=[P_TPM.time, P_TPM.lev, lat_gc_AK, lon_gc_AK], dims=['time','lev' ,'lat','lon'])

# remove 0
BXHEIGHT_TPM_empty[BXHEIGHT_TPM_empty==0]=np.nan
# pack in xarray
BXHEIGHT_TPM = xr.DataArray(BXHEIGHT_TPM_empty, coords=[P_TPM.time, P_TPM.lev, lat_gc_AK, lon_gc_AK], dims=['time','lev' ,'lat','lon'])

# remove 0
gc_ND_HCHO_TPMgrids_empty[gc_ND_HCHO_TPMgrids_empty==0]=np.nan
# pack in xarray
gc_ND_HCHO_TPMgrids = xr.DataArray(gc_ND_HCHO_TPMgrids_empty, coords=[P_TPM.time, P_TPM.lev, lat_gc_AK, lon_gc_AK], dims=['time','lev' ,'lat','lon'])








# TPM Number Density in each TROPOMI levels
NumberDensity_TPM = P_TPM*6.02e23/(8.314*T_TPM)*1e-6*apriori_TPM

# TROPOMI AMF
tpm_BXN = NumberDensity_TPM*BXHEIGHT_TPM
sum_tpm_BXN = tpm_BXN.sum(dim='lev')
AMF_TPM = (tpm_BXN /sum_tpm_BXN * AveragingKernel_TPM * AMF_HCHO_trop_TPM ).sum(dim='lev')
# FillValue -> NaN
AMF_TPM.values[AMF_TPM==0]=np.nan




# GEOS-Chem AMF
gc_BXN = gc_ND_HCHO_TPMgrids*BXHEIGHT_TPM
sum_gc_BXN = gc_BXN.sum(dim='lev')
AMF_GC = (gc_BXN /sum_gc_BXN * AveragingKernel_TPM * AMF_HCHO_trop_TPM ).sum(dim='lev')
# FillValue -> NaN
AMF_GC.values[AMF_GC==0]=np.nan



# write GCAMF
AMF_GC.to_netcdf(runname_dir+"AMF_GC_"+YYYY[2:4]+".revised.2.nc", engine="netcdf4")
# write TPM AMF
AMF_TPM.to_netcdf(runname_dir+"AMF_TPMcalc_"+YYYY[2:4]+".revised.2.nc", engine="netcdf4")
# write Vertical profile of GC in TPM grid
gc_BXN.to_netcdf(runname_dir+"HCHO_apriori_GC_"+YYYY[2:4]+".revised.2.nc", engine="netcdf4")



