##################################################################
#                                                                #
#                              MERRA2                            #
#                                                                #
################################################################## 
import xarray as xr
import numpy as np

merra19 = xr.open_mfdataset('../../../../storage/edb/MERRA2/19*/*.nc')
merra20 = xr.open_mfdataset('../../../../storage/edb/MERRA2/20*/*.nc')
merra = xr.concat([merra19, merra20], dim='time')
merra = merra.sel(time=slice('1980','2018'))
merra = merra.where(merra.lat >= 17, drop=True)
merra = merra.where(merra.lat <= 33, drop=True)
merra = merra.where(merra.lon >= -100, drop=True)
merra = merra.where(merra.lon <= -80, drop=True)
merra['WS_10'] = np.sqrt(merra.U10M**2+merra.V10M**2)

#Interpolacion a 90m
merra['WS_90'] = merra['WS_10']*((90/10)**0.11)

#Interpolacion a 119m
merra['WS_119'] = merra['WS_10']*((119/10)**0.11)

#Promedio global
merra_promedio_global = merra.mean(dim='time')

merra_global_90 = np.asarray(merra_promedio_global.WS_90)
np.savetxt('merra_global_90.csv', merra_global_90)

merra_global_119 = np.asarray(merra_promedio_global.WS_119)
np.savetxt('merra_global_119.csv', merra_global_119)

#Promedios anuales
merra_anual = merra.resample(time='Y').mean()

years = np.arange(0,39,1)

merra_lat = np.asarray(merra_anual.lat)
np.savetxt('merra_lat.csv', merra_lat)

merra_lon = np.asarray(merra_anual.lon)
np.savetxt('merra_lon.csv', merra_lon)

###### Datos a 90M
merra_ws_90 = np.asarray(merra.WS_90)

#Save csv
for i in years:
	np.savetxt('merra_anual_mean_90m_{}.csv'.format(1980+i), merra_ws_90[i])


##### Datos a 119M
#merra_ws_119 = np.asarray(merra_anual.WS_119)
merra_ws_119 = np.asarray(merra.WS_119)

#Save csv
for i in years:
	np.savetxt('merra_anual_mean_119m_{}.csv'.format(1980+i), merra_ws_119[i])


###################################################################
#                                                                 #
#                             ERA5                                #
#                                                                 #
###################################################################
import xarray as xr
import numpy as np

era = xr.open_mfdataset('../../../../storage/edb/Era-5/*.nc')
era = era.sel(time=slice('1980','2018'))
era = era.where(era.latitude >= 17, drop=True)
era = era.where(era.latitude <= 33, drop=True)
era = era.where(era.longitude >= 260, drop=True)
era = era.where(era.longitude <= 280, drop=True)
era['WS_10'] = np.sqrt(era.u10**2+era.v10**2)

#Interpolacion a 90m
era['WS_90'] = era['WS_10']*((90/10)**0.11)

#Interpolacion a 119m
era['WS_119'] = era['WS_10']*((119/10)**0.11)

#Promedio global
era_promedio_global = era.mean(dim='time')

era_global_90 = np.asarray(era_promedio_global.WS_90)
np.savetxt('era_global_90.csv', era_global_90)


era_global_119 = np.asarray(era_promedio_global.WS_119)
np.savetxt('era_global_119.csv', era_global_119)

era_anual = era.resample(time='Y').mean()

years = np.arange(0,39,1)

era_lat = np.asarray(era_anual.latitude)
np.savetxt('era_lat.csv', era_lat)

era_lon = np.asarray(era_anual.longitude)
np.savetxt('era_lon.csv', era_lon)

###### Datos a 90m
#era_ws_90 = np.asarray(era_anual.WS_90)
era_ws_90 = np.asarray(era.WS_90)

#Save CSV
for i in years:
	np.savetxt('era_anual_mean_90m_{}.csv'.format(1980+i), era_ws_90[i])

###### Datos a 119m
#era_ws_119 = np.asarray(era_anual.WS_119)
era_ws_119 = np.asarray(era.WS_119)

#Save CSV
for i in years:
	np.savetxt('era_anual_mean_119m_{}.csv'.format(1980+i), era_ws_119[i])


############################################################################
#                                                                          #
#                                  Curva DTU                               #
#                                                                          #
############################################################################

#Definir curva
velocidades = np.linspace(0,28,29)
power = [0, 0, 0, 0, 280.2, 799.1, 1532.7, 2506.1, 3730.7, 5311.8, 7286.5, 9698.3, 10639.1, 10648.5, 10639.3, 10683.7, 10642.0, 10640.0, 10639.9, 10652.8, 10646.2, 10644.0, 10641.2, 10639.5, 10643.6, 10635.7, 0, 0, 0]

fpt1 = np.poly1d(np.polyfit(np.linspace(2, 12, 11),[0, 0, 280.2, 799.1, 1532.7, 2506.1, 3730.7, 5311.8, 7286.5, 9698.3, 10639.1], deg = 6))

def curva_aero(data):
    vel_turb = []
    for i in data:
        if i < 4:
            vel_turb.append(0)
        elif i >= 4 and i < 11.7:
            vel_turb.append(fpt1(i))
        elif i >=11.7 and i <= 25:
            vel_turb.append(10640)
        else:
            vel_turb.append(0)
    return(vel_turb)

#Evaluar curva
potencia = []
for i in range(65):
    for j in range(81):
        potencia.append(sum(curva_aero(era_ws_119[:,i,j])))
        print(i,j)


potencia = np.asarray(potencia)
potencia = potencia.reshape(65,81)
np.savetxt('potencia_era_dtu.csv', potencia)



############################################################################
#                                                                          #
#                                  Curva NREL                              #
#                                                                          #
############################################################################

velocidades = np.linspace(0,28,29)

power_nrel = [0,0,0,40.5,177.7,403.9,737.6,1187.2,1771.1,2518.6,3448.4,4562.5,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,0,0,0]

fpt2 = np.poly1d(np.polyfit(np.arange(1,13,1),[0,0,40.5,177.7,403.9,737.6,1187.2,1771.1,2518.6,3448.4,4562.5,5000], deg=6))

def curva_aero_nrel(data):
    vel_turb = []
    for i in data:
        if i < 0:
            vel_turb.append(0)
        elif i >= 3 and i <11.7:
            vel_turb.append(fpt2(i))
        elif i >= 11.7 and i <= 25:
            vel_turb.append(5000)
        else:
            vel_turb.append(0)
    return(vel_turb)

#Evaluar curva
potencia = []
for i in range(33):
    for j in range(33):
        potencia.append(sum(curva_aero_nrel(merra_ws_90[:,i,j])))
        print(i,j)


potencia = np.asarray(potencia)
potencia = potencia.reshape(33,33)
np.savetxt('potencia_merra_nrel.csv', potencia)



