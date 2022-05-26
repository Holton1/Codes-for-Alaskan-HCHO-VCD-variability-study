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
import glob
import cmaps
import sys

def file_list(dirname, ext='.csv'):
    import os
    return list(filter(
        lambda filename: os.path.splitext(filename)[1] == ext,
        os.listdir(dirname)))

import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

mpl.rcParams['font.sans-serif'] = "Arial"
mpl.rcParams['font.family'] = "sans-serif"
plt.rcParams['figure.figsize'] = (10, 6)
mpl.rc('xtick', labelsize=18) 
mpl.rc('ytick', labelsize=18)
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.titlesize'] = 18

import rum

# new colormap
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
cmap = cmaps.GHRSST_anomaly
newcolors=cmap(np.linspace(0, 1, 256))
GHRSST_pos = ListedColormap(newcolors[128:])


def proplot_mjja(var,tmp_title,vminmax,unitname,cmapname,cbposition):

    ############################    show   GC   restults    ################################
    import numpy as np
    import Ngl,Nio
    import matplotlib.pyplot as plt
    from matplotlib.dates import DateFormatter
    import matplotlib as mpl

    import cartopy.feature as cfeature
    import cartopy.crs as ccrs
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    import matplotlib.cm as cm
    monname = ["May","June","July","August"]
    str_ipanelindex = ["a","b","c",'d','e','f','g','h','i','j','k','l']
    
    varnum = len(var)
    nrow = varnum
    ncol = 4
    
    panelwidth = 3.5
    figsize = (panelwidth*ncol , panelwidth*nrow)
    fig = plt.figure(figsize=figsize)

    ipanel = 0
    ipanelindex = np.arange(1,1+nrow*ncol,1)
    
    
    for ivar in range(0,varnum):

        # monthly mean variable
        mm_HCHOVCD = var[ivar].resample(time="1MS").mean(dim="time")   # 5678
        iseason = 0
        
        for imon in range(5,9): 

            vmin=vminmax[ivar][0]
            vmax=vminmax[ivar][1]
            cmap = cmapname[ivar]
            norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

            #####################################################################################
            ###################  GEOSCHEM HCHO VCD regridded noon (time mean)
            #####################################################################################

            ax = fig.add_subplot(nrow,ncol,ipanelindex[ipanel], projection=ccrs.PlateCarree(), aspect='auto')

            # Select some data to plot
            im = ax.pcolormesh(mm_HCHOVCD.lon, mm_HCHOVCD.lat, mm_HCHOVCD[imon-5,:,:] ,cmap=cmap,norm=norm)

            # add colorbar
            mappable = cm.ScalarMappable(cmap=cmap)
            mappable.set_array([])
            mappable.set_clim(vmin=vmin, vmax=vmax)
                    
            # add lat lon label
#             ax.set_xticks([-160,-150,-140,-130], crs=ccrs.PlateCarree())
#             if ipanelindex[ipanel]==1:
#                 ax.set_yticks([50,55,60,65,70,75], crs=ccrs.PlateCarree())
#             lon_formatter = LongitudeFormatter(zero_direction_label=True)
#             lat_formatter = LatitudeFormatter()
#             ax.xaxis.set_major_formatter(lon_formatter)
#             ax.yaxis.set_major_formatter(lat_formatter)


#             ax.annotate( str_ipanelindex[ipanel], xy=(0.07, 0.08), xycoords="axes fraction", fontweight='bold', fontsize = 30 )
            if ivar==0:
                plt.title(monname[imon-5], fontsize=20)
#             else:
#                 ax.figtext(0.1,0.9,tmp_title[ivar])
            if iseason==3:  
                ax.annotate(tmp_title[ivar], xy=(0.03-3, 0.9), xycoords="axes fraction",fontsize=20, fontweight='bold')

#             plt.title(monname[imon-5]+" "+tmp_title[ivar])    
            ax.add_feature(cfeature.COASTLINE)
            ax.add_feature(cfeature.BORDERS)        
            land_50m = cfeature.NaturalEarthFeature('physical', 'coastline', '110m',
                                                    edgecolor='k',
                                                    facecolor='none')
            ax.add_feature(land_50m)
            
            ipanel+=1
            iseason+=1

        if cbposition == "none":
            continue
        else:
            import matplotlib.ticker                         # here's where the formatter is
            cbformat = matplotlib.ticker.ScalarFormatter()   # create the formatter
            cbformat.set_powerlimits((-2,2))                 # set the limits for sci. not.
    #         position=fig.add_axes([0.91, 0.65, 0.005, 0.23])  # original good scale
            position=fig.add_axes( cbposition )
            cb = fig.colorbar(im,cax=position,orientation='vertical',format=cbformat)
            cb.set_label(unitname[ivar], fontsize=20)


        fig.tight_layout()#调整整体空白 
        plt.subplots_adjust(wspace =0.01, hspace =0.01)#调整子图间距
        # plt.savefig(outdir+'???',dpi=300)
        


runname = sys.argv[1]
YYYY = sys.argv[2]

jndirGC="/import/GREENING/tzhao/jndata/GEOS-Chem/"
jndirGC_WF = jndirGC+"MEGANon_PFTol_WF"+"/"
jndirGC_runname = jndirGC+runname+"/"









# read in raw Emission data

filetype = "HEMCO_diagnostics"


ds_EmisXYLE_BioBurn = []
ds_EmisTOLU_BioBurn = []
ds_EmisSOAP_BioBurn = []
ds_EmisPRPE_BioBurn = []
ds_EmisOCPO_BioBurn = []
ds_EmisOCPI_BioBurn = []
ds_EmisMTPA_BioBurn = []
ds_EmisMOH_BioBurn = []
ds_EmisMEK_BioBurn = []
ds_EmisEOH_BioBurn = []
ds_EmisCO_BioBurn = []
ds_EmisCH4_BioBurn = []
ds_EmisCH2O_BioBurn = []
ds_EmisC3H8_BioBurn = []
ds_EmisC2H6_BioBurn = []
ds_EmisBENZ_BioBurn = []
ds_EmisBCPO_BioBurn = []
ds_EmisBCPI_BioBurn = []
ds_EmisALK4_BioBurn = []
ds_EmisALD2_BioBurn = []
ds_EmisACET_BioBurn = []  



for MM in range(5,9):
    
    dir_GCout = "/import/GREENING/tzhao/"+str(runname)+"/"+str(YYYY)+"/nest"+str(MM).zfill(2)+"_merra2_2x25_tropchem/OutputDir/"    
    filelist = list(map(os.path.basename, sorted( glob.glob(dir_GCout + filetype + "*" + str(YYYY)+str(MM).zfill(2) +"*")) ) )
    
    for fname in filelist:
        
        ds = xr.open_dataset( dir_GCout + fname , engine="netcdf4")
        
        ds_EmisXYLE_BioBurn.append(ds["EmisXYLE_BioBurn"])
        ds_EmisTOLU_BioBurn.append(ds["EmisTOLU_BioBurn"])
        ds_EmisSOAP_BioBurn.append(ds["EmisSOAP_BioBurn"])
        ds_EmisPRPE_BioBurn.append(ds["EmisPRPE_BioBurn"])
        ds_EmisOCPO_BioBurn.append(ds["EmisOCPO_BioBurn"])
        ds_EmisOCPI_BioBurn.append(ds["EmisOCPI_BioBurn"])
        ds_EmisMTPA_BioBurn.append(ds["EmisMTPA_BioBurn"])
        ds_EmisMOH_BioBurn.append(ds["EmisMOH_BioBurn"])
        ds_EmisMEK_BioBurn.append(ds["EmisMEK_BioBurn"])
        ds_EmisEOH_BioBurn.append(ds["EmisEOH_BioBurn"])
        ds_EmisCO_BioBurn.append(ds["EmisCO_BioBurn"])
        ds_EmisCH4_BioBurn.append(ds["EmisCH4_BioBurn"])
        ds_EmisCH2O_BioBurn.append(ds["EmisCH2O_BioBurn"])
        ds_EmisC3H8_BioBurn.append(ds["EmisC3H8_BioBurn"])
        ds_EmisC2H6_BioBurn.append(ds["EmisC2H6_BioBurn"])
        ds_EmisBENZ_BioBurn.append(ds["EmisBENZ_BioBurn"])
        ds_EmisBCPO_BioBurn.append(ds["EmisBCPO_BioBurn"])
        ds_EmisBCPI_BioBurn.append(ds["EmisBCPI_BioBurn"])
        ds_EmisALK4_BioBurn.append(ds["EmisALK4_BioBurn"])
        ds_EmisALD2_BioBurn.append(ds["EmisALD2_BioBurn"])
        ds_EmisACET_BioBurn.append(ds["EmisACET_BioBurn"])

        
    
    

EmisXYLE_BioBurn = xr.concat(ds_EmisXYLE_BioBurn, 'time')
EmisTOLU_BioBurn = xr.concat(ds_EmisTOLU_BioBurn, 'time')
EmisSOAP_BioBurn = xr.concat(ds_EmisSOAP_BioBurn, 'time')
EmisPRPE_BioBurn = xr.concat(ds_EmisPRPE_BioBurn, 'time')
EmisOCPO_BioBurn = xr.concat(ds_EmisOCPO_BioBurn, 'time')
EmisOCPI_BioBurn = xr.concat(ds_EmisOCPI_BioBurn, 'time')
EmisMTPA_BioBurn = xr.concat(ds_EmisMTPA_BioBurn, 'time')
EmisMOH_BioBurn = xr.concat(ds_EmisMOH_BioBurn, 'time')
EmisMEK_BioBurn = xr.concat(ds_EmisMEK_BioBurn, 'time')
EmisEOH_BioBurn = xr.concat(ds_EmisEOH_BioBurn, 'time')
EmisCO_BioBurn = xr.concat(ds_EmisCO_BioBurn, 'time')
EmisCH4_BioBurn = xr.concat(ds_EmisCH4_BioBurn, 'time')
EmisCH2O_BioBurn = xr.concat(ds_EmisCH2O_BioBurn, 'time')
EmisC3H8_BioBurn = xr.concat(ds_EmisC3H8_BioBurn, 'time')
EmisC2H6_BioBurn = xr.concat(ds_EmisC2H6_BioBurn, 'time')
EmisBENZ_BioBurn = xr.concat(ds_EmisBENZ_BioBurn, 'time')
EmisBCPO_BioBurn = xr.concat(ds_EmisBCPO_BioBurn, 'time')
EmisBCPI_BioBurn = xr.concat(ds_EmisBCPI_BioBurn, 'time')
EmisALK4_BioBurn = xr.concat(ds_EmisALK4_BioBurn, 'time')
EmisALD2_BioBurn = xr.concat(ds_EmisALD2_BioBurn, 'time')
EmisACET_BioBurn = xr.concat(ds_EmisACET_BioBurn, 'time')



data_bioburn = xr.Dataset( 
    {
        'EmisXYLE_BioBurn': ( ('time', 'lat',  'lon'),  EmisXYLE_BioBurn[:,5:-5,4:-4].values ),\
        'EmisTOLU_BioBurn': ( ('time', 'lat',  'lon'),  EmisTOLU_BioBurn[:,5:-5,4:-4].values ),\
        'EmisSOAP_BioBurn': ( ('time', 'lat',  'lon'),  EmisSOAP_BioBurn[:,5:-5,4:-4].values ),\
        'EmisPRPE_BioBurn': ( ('time', 'lat',  'lon'),  EmisPRPE_BioBurn[:,5:-5,4:-4].values ),\
        'EmisOCPO_BioBurn': ( ('time', 'lat',  'lon'),  EmisOCPO_BioBurn[:,5:-5,4:-4].values ),\
        'EmisOCPI_BioBurn': ( ('time', 'lat',  'lon'),  EmisOCPI_BioBurn[:,5:-5,4:-4].values ),\
        'EmisMTPA_BioBurn': ( ('time', 'lat',  'lon'),  EmisMTPA_BioBurn[:,5:-5,4:-4].values ),\
        'EmisMOH_BioBurn': ( ('time', 'lat',  'lon'),  EmisMOH_BioBurn[:,5:-5,4:-4].values ),\
        'EmisMEK_BioBurn': ( ('time', 'lat',  'lon'),  EmisMEK_BioBurn[:,5:-5,4:-4].values ),\
        'EmisEOH_BioBurn': ( ('time', 'lat',  'lon'),  EmisEOH_BioBurn[:,5:-5,4:-4].values ),\
        'EmisCO_BioBurn': ( ('time', 'lat',  'lon'),  EmisCO_BioBurn[:,5:-5,4:-4].values ),\
        'EmisCH4_BioBurn': ( ('time', 'lat',  'lon'),  EmisCH4_BioBurn[:,5:-5,4:-4].values ),\
        'EmisCH2O_BioBurn': ( ('time', 'lat',  'lon'),  EmisCH2O_BioBurn[:,5:-5,4:-4].values ),\
        'EmisC3H8_BioBurn': ( ('time', 'lat',  'lon'),  EmisC3H8_BioBurn[:,5:-5,4:-4].values ),\
        'EmisC2H6_BioBurn': ( ('time', 'lat',  'lon'),  EmisC2H6_BioBurn[:,5:-5,4:-4].values ),\
        'EmisBENZ_BioBurn': ( ('time', 'lat',  'lon'),  EmisBENZ_BioBurn[:,5:-5,4:-4].values ),\
        'EmisBCPO_BioBurn': ( ('time', 'lat',  'lon'),  EmisBCPO_BioBurn[:,5:-5,4:-4].values ),\
        'EmisBCPI_BioBurn': ( ('time', 'lat',  'lon'),  EmisBCPI_BioBurn[:,5:-5,4:-4].values ),\
        'EmisALK4_BioBurn': ( ('time', 'lat',  'lon'),  EmisALK4_BioBurn[:,5:-5,4:-4].values ),\
        'EmisALD2_BioBurn': ( ('time', 'lat',  'lon'),  EmisALD2_BioBurn[:,5:-5,4:-4].values ),\
        'EmisACET_BioBurn': ( ('time', 'lat',  'lon'),  EmisACET_BioBurn[:,5:-5,4:-4].values )
    },\
    
    coords = 
    {
        'time': EmisXYLE_BioBurn[:,5:-5,4:-4].time,\
        'lat': EmisXYLE_BioBurn[:,5:-5,4:-4].lat,\
        'lon': EmisXYLE_BioBurn[:,5:-5,4:-4].lon   
    }
)


# write emission
data_bioburn.to_netcdf(jndirGC_runname+"Emis_BioBurn_"+YYYY[-2:]+".nc", engine="netcdf4")

print("write Emis_BioBurn done")




