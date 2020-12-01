import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
import pandas as pd

velocidades = np.linspace(0,28,29)
power = [0, 0, 0, 0, 280.2, 799.1, 1532.7, 2506.1, 3730.7, 5311.8, 7286.5, 9698.3, 10639.1, 
         10648.5,  10639.3, 10683.7, 10642.0, 10640.0, 10639.9, 10652.8, 10646.2, 10644.0, 
         10641.2, 10639.5, 10643.6, 10635.7, 0, 0, 0]

turbine = pd.DataFrame()
turbine['Velocidades'] = velocidades
turbine['Power'] = power

fpt1 = np.poly1d(np.polyfit(np.linspace(2, 12, 11), 
                 [0, 0, 280.2, 799.1, 1532.7, 2506.1, 3730.7, 5311.8, 7286.5, 9698.3, 10639.1], deg = 6))

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

merra_anual = xr.open_mfdataset('./merra_2018/*.nc')
merra_anual = merra_anual.where(merra_anual.lat >= 17, drop = True)
merra_anual = merra_anual.where(merra_anual.lon >= -99.375, drop = True)

merra_anual['WS50'] = (merra_anual['U50M']**2 + merra_anual['V50M']**2)**(1/2)
merra_anual['WSrot'] = merra_anual['WS50']*(119/50)**(1/7)

time = merra_anual.variables['time']
lons = merra_anual.variables['lon']
lats = merra_anual.variables['lat']

groups = dict(list(merra_anual.groupby('time.season')))

WSrot_anual = np.asarray(merra_anual.variables['WSrot'])
WSrot_inv  =  np.asarray(groups['DJF'].variables['WSrot'])
WSrot_prim =  np.asarray(groups['MAM'].variables['WSrot'])
WSrot_vera =  np.asarray(groups['JJA'].variables['WSrot'])
WSrot_oto  =  np.asarray(groups['SON'].variables['WSrot'])

print('empieza anual')
potencia_anual = []
for i in range(33):
    for j in range(32):
        potencia_anual.append(sum(curva_aero(WSrot_anual[:,i,j])))
        
potencia_anual = np.asarray(potencia_anual)
potencia_anual = potencia_anual.reshape(33,32)

print('empieza invierno')
potencia_inv = []
for i in range(33):
    for j in range(32):
        potencia_inv.append(sum(curva_aero(WSrot_inv[:,i,j])))
        
potencia_inv = np.asarray(potencia_inv)
potencia_inv = potencia_inv.reshape(33,32)

print('empieza primavera')
potencia_prim = []
for i in range(33):
    for j in range(32):
        potencia_prim.append(sum(curva_aero(WSrot_prim[:,i,j])))
        
potencia_prim = np.asarray(potencia_prim)
potencia_prim = potencia_prim.reshape(33,32)

print('empieza verano')
potencia_vera = []
for i in range(33):
    for j in range(32):
        potencia_vera.append(sum(curva_aero(WSrot_vera[:,i,j])))
        
potencia_vera = np.asarray(potencia_vera)
potencia_vera = potencia_vera.reshape(33,32)

print('empieza otoño')
potencia_oto = []
for i in range(33):
    for j in range(32):
        potencia_oto.append(sum(curva_aero(WSrot_oto[:,i,j])))
        
potencia_oto = np.asarray(potencia_oto)
potencia_oto = potencia_oto.reshape(33,32)


###########################################
###########################################

## gebco = xr.open_dataset('./gebco/GEBCO_2019_-120.3189_34.6763_-79.658_15.1747_mexico.nc')
## gebco = gebco.where(gebco.lat >= 17, drop = True)
## gebco = gebco.where(gebco.lon >= -99.375, drop = True)

## lonsg = gebco.variables['lon']
## latsg = gebco.variables['lat']
## long, latg = np.meshgrid(lonsg, latsg)

#elev = gebco.elevation.values
## elevuns = np.asarray(gebco.variables['elevation'])

###########################################
###########################################
print('empieza grafica')
plt.figure(figsize=(8,16))
grid = plt.GridSpec(4, 2, hspace=0.3)

plt.subplot(grid[0:2,0:2])
map=Basemap(projection='merc',llcrnrlon=-99.375,llcrnrlat=17,urcrnrlon=-80.0,urcrnrlat=33,resolution='i')
# projection, lat/lon extents and resolution of polygons to draw
# resolutions: c - crude, l - low, i - intermediate, h - high, f - full
lon, lat = np.meshgrid(lons, lats)
xi, yi = map(lon, lat)

# plot wind speed and direction
cs = map.contourf(xi,yi,(potencia_anual)/(10640*8760)*100, [0,10,20,30,40,50,60,70], cmap=cm.rainbow)
#cs.set_edgecolor('face')
#map.quiver(xi, yi, U2M_nans[0,:,:], V2M_nans[0,:,:], scale=600, color='k')

# Add Grid Lines
#map.drawparallels(np.linspace(17,33,33), labels = [1,0,0,0],fontsize=10, linewidth=0.3)
#map.drawmeridians(np.linspace(-99.375,80,38), labels=[0,0,0,1], fontsize=10, linewidth=0.3)

# Add Coastlines, States, and Country Boundaries
map.drawcoastlines(color='k', linewidth=0.5)
map.drawstates(color='k', linewidth=0.5)
map.drawcountries(color='k', linewidth=0.5)
#map.fillcontinents(color = 'coral')
#map.drawstates(color='grey', linewidth=0.5)
#map.drawcountries(color='k', linewidth=0.5)
#map.arcgisimage(service = 'World_Shaded_Relief', xpixels = 2000)

# Add Colorbar
cbar = map.colorbar(cs, location='bottom', pad="5%")
cbar.set_label('Capacity Factor [%]')
cbar.ax.tick_params(labelsize=8)
plt.title('Anual 2018', fontsize = 12)

###########################################################
#################################
plt.subplot(grid[2,0])
map = Basemap(projection='merc',llcrnrlon=-99.375,llcrnrlat=17,urcrnrlon=-80.0,urcrnrlat=33,resolution='i')
# projection, lat/lon extents and resolution of polygons to draw
# resolutions: c - crude, l - low, i - intermediate, h - high, f - full
lon, lat = np.meshgrid(lons, lats)
xi, yi = map(lon, lat)

# plot wind speed and direction
cs = map.contourf(xi,yi,(potencia_inv)/(10640*2160)*100, [0,10,20,30,40,50,60,70], cmap=cm.rainbow)
#cs.set_edgecolor('face')
#map.quiver(xi, yi, U2M_nans[0,:,:], V2M_nans[0,:,:], scale=600, color='k')

# Add Grid Lines
#map.drawparallels(np.linspace(17,33,33), labels = [1,0,0,0],fontsize=10, linewidth=0.3)
#map.drawmeridians(np.linspace(-99.375,80,38), labels=[0,0,0,1], fontsize=10, linewidth=0.3)

# Add Coastlines, States, and Country Boundaries
map.drawcoastlines(color='k', linewidth=0.5)
map.drawstates(color='k', linewidth=0.5)
map.drawcountries(color='k', linewidth=0.5)
#map.fillcontinents(color = 'coral')
#map.drawstates(color='grey', linewidth=0.5)
#map.drawcountries(color='k', linewidth=0.5)
#map.arcgisimage(service = 'World_Shaded_Relief', xpixels = 2000)

# Add Colorbar
cbar = map.colorbar(cs, location='bottom', pad="5%")
cbar.set_label('Capacity Factor [%]')
cbar.ax.tick_params(labelsize=8)
plt.title('Winter', fontsize = 12)

###########################################################

plt.subplot(grid[2,1])
map = Basemap(projection='merc',llcrnrlon=-99.375,llcrnrlat=17,urcrnrlon=-80.0,urcrnrlat=33,resolution='i')
# projection, lat/lon extents and resolution of polygons to draw
# resolutions: c - crude, l - low, i - intermediate, h - high, f - full
lon, lat = np.meshgrid(lons, lats)
xi, yi = map(lon, lat)

# plot wind speed and direction
cs = map.contourf(xi,yi,(potencia_prim)/(10640*2208)*100, [0,10,20,30,40,50,60,70], cmap=cm.rainbow)
#cs.set_edgecolor('face')
#map.quiver(xi, yi, U2M_nans[0,:,:], V2M_nans[0,:,:], scale=600, color='k')

# Add Grid Lines
#map.drawparallels(np.linspace(17,33,33), labels = [1,0,0,0],fontsize=10, linewidth=0.3)
#map.drawmeridians(np.linspace(-99.375,80,38), labels=[0,0,0,1], fontsize=10, linewidth=0.3)

# Add Coastlines, States, and Country Boundaries
map.drawcoastlines(color='k', linewidth=0.5)
map.drawstates(color='k', linewidth=0.5)
map.drawcountries(color='k', linewidth=0.5)
#map.fillcontinents(color = 'coral')
#map.drawstates(color='grey', linewidth=0.5)
#map.drawcountries(color='k', linewidth=0.5)
#map.arcgisimage(service = 'World_Shaded_Relief', xpixels = 2000)

# Add Colorbar
cbar = map.colorbar(cs, location='bottom', pad="5%")
cbar.set_label('Capacity Factor [%]')
cbar.ax.tick_params(labelsize=8)
plt.title('Spring')

########################################

plt.subplot(grid[3,0])
map = Basemap(projection='merc',llcrnrlon=-99.375,llcrnrlat=17,urcrnrlon=-80.0,urcrnrlat=33,resolution='i')
# projection, lat/lon extents and resolution of polygons to draw
# resolutions: c - crude, l - low, i - intermediate, h - high, f - full
lon, lat = np.meshgrid(lons, lats)
xi, yi = map(lon, lat)

# plot wind speed and direction
cs = map.contourf(xi,yi,(potencia_vera)/(10640*2184)*100, [0,10,20,30,40,50,60,70], cmap=cm.rainbow)
#cs.set_edgecolor('face')
#map.quiver(xi, yi, U2M_nans[0,:,:], V2M_nans[0,:,:], scale=600, color='k')

# Add Grid Lines
#map.drawparallels(np.linspace(17,33,33), labels = [1,0,0,0],fontsize=10, linewidth=0.3)
#map.drawmeridians(np.linspace(-99.375,80,38), labels=[0,0,0,1], fontsize=10, linewidth=0.3)

# Add Coastlines, States, and Country Boundaries
map.drawcoastlines(color='k', linewidth=0.5)
map.drawstates(color='k', linewidth=0.5)
map.drawcountries(color='k', linewidth=0.5)
#map.fillcontinents(color = 'coral')
#map.drawstates(color='grey', linewidth=0.5)
#map.drawcountries(color='k', linewidth=0.5)
#map.arcgisimage(service = 'World_Shaded_Relief', xpixels = 2000)

# Add Colorbar
cbar = map.colorbar(cs, location='bottom', pad="5%")
cbar.set_label('Capacity Factor [%]')
cbar.ax.tick_params(labelsize=8)
plt.title('Summer')

##############################################

plt.subplot(grid[3,1])
map = Basemap(projection='merc',llcrnrlon=-99.375,llcrnrlat=17,urcrnrlon=-80.0,urcrnrlat=33,resolution='i')
# projection, lat/lon extents and resolution of polygons to draw
# resolutions: c - crude, l - low, i - intermediate, h - high, f - full
lon, lat = np.meshgrid(lons, lats)
xi, yi = map(lon, lat)

# plot wind speed and direction
cs = map.contourf(xi,yi,(potencia_oto)/(10640*2208)*100, [0,10,20,30,40,50,60,70], cmap=cm.rainbow)
#cs.set_edgecolor('face')
#map.quiver(xi, yi, U2M_nans[0,:,:], V2M_nans[0,:,:], scale=600, color='k')

# Add Grid Lines
#map.drawparallels(np.linspace(17,33,33), labels = [1,0,0,0],fontsize=10, linewidth=0.3)
#map.drawmeridians(np.linspace(-99.375,80,38), labels=[0,0,0,1], fontsize=10, linewidth=0.3)

# Add Coastlines, States, and Country Boundaries
map.drawcoastlines(color='k', linewidth=0.5)
map.drawstates(color='k', linewidth=0.5)
map.drawcountries(color='k', linewidth=0.5)
#map.fillcontinents(color = 'coral')
#map.drawstates(color='grey', linewidth=0.5)
#map.drawcountries(color='k', linewidth=0.5)
#map.arcgisimage(service = 'World_Shaded_Relief', xpixels = 2000)



# Add Colorbar
cbar = map.colorbar(cs, location='bottom', pad="5%")
cbar.set_label('Capacity Factor [%]')
cbar.ax.tick_params(labelsize=8)
plt.title('Autumn')

#plt.savefig('capacity_factor18_nofill_rel.pdf', dpi=360)

