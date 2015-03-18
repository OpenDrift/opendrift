import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


members=11
times=100
# Define a constant grid:
latmin = 57.0
lonmin = 8.5
latmax = 59.0
lonmax = 11.0
print latmin, lonmin, latmax, lonmax
 
dlat = 0.01 #these can be found automatically from the lat/lon vaules themselves
dlon = 0.01
 
lon2d, lat2d = np.meshgrid(np.linspace(lonmin, lonmax, int((lonmax-lonmin)/dlon)+1), 
                           np.linspace(latmin, latmax, int((latmax-latmin)/dlat)+1))


def read_member(t,member):
    txt = pd.read_csv('eps_out/tstep_'+str(t)+'_member'+str(member)+'.dat').values
    lat         = np.empty(len(txt),dtype=np.float)
    lon         = np.empty(len(txt),dtype=np.float)
    id          = np.empty(len(txt),dtype=np.float)
    for j in range(len(id)):
        id[j]  = txt[j][0]
        lat[j] = txt[j][1]
        lon[j] = txt[j][2]
    return id, lat, lon


for t in range(times):
    id =[]
    lat = []
    lon =[]
    for n in range(members):
        a, b, c = read_member(t,n)
        id.append(a)
        lat.append(b)
        lon.append(c)

    #Mask nans
    lat = np.ma.masked_invalid(lat)
    lon = np.ma.masked_invalid(lon)
    
#     latmin = np.min(lat)
#     lonmin = np.min(lon)
#     latmax = np.max(lat)
#     lonmax = np.max(lon)
#     print latmin, lonmin, latmax, lonmax
#      
#     dlat = 0.01 #these can be found automatically from the lat/lon vaules themselves
#     dlon = 0.01
#      
#     lon2d, lat2d = np.meshgrid(np.linspace(lonmin, lonmax, int((lonmax-lonmin)/dlon)+1), 
#                                np.linspace(latmin, latmax, int((latmax-latmin)/dlat)+1))
    # lon2d, lat2d = np.meshgrid(np.array(range(lonmin, lonmax, dlon)), 
    #                            np.array(range(latmin, latmax, dlat)))
    count = np.zeros_like(lon2d)
    
    xp = np.hstack(lon)
    yp = np.hstack(lat)
    for i,j in map(None,xp,yp):
        y = (np.abs(lat2d[:,0]-j)).argmin()
        x = (np.abs(lon2d[0,:]-i)).argmin()
    #     y = np.where(lat2d[:,0]==j)[0]
    #     x = np.where(lon2d[0,:]==i)[0]
        #print x, y
        count[y,x]+=1
    #print count
    
    count[0,0]=0 #HACK! must deal with the NaNs in a better way in the future...
    
    plt.clf()
    #plt.plot(lon,lat,'r.')
    ###plt.contourf(lon2d,lat2d,(count/np.max(count))*100,np.arange(5,100,1),cmap='PuRd', extend='max') 
    #plt.pcolormesh(lon2d,lat2d,(count/np.max(count))*100,vmin=5,vmax=100,cmap='PuRd') 
    plt.pcolormesh(lon2d,lat2d,count,vmin=5,vmax=10,cmap='PuRd') 
    #plt.contour(lon2d,lat2d,count,[1,10,100]) 
    # plt.pcolormesh(lon2d,lat2d,count, 
    #                norm=LogNorm(vmin=0, vmax=count.max()), cmap='PuBu_r')
    plt.colorbar()
    name = 'eps_out/time%02i.png' % t
    plt.savefig(name,dpi=150,facecolor='white')
    #plt.show()
    print 'done '  +str(t)

print 'done'