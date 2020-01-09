import netCDF4 as NC
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

dat = NC.Dataset("project3.nc","r",format="NETCDF4")

#print(dat.__dict__)
#print(dat.variables['rho_charge_density'])
x = dat.variables['x']
h0 = dat.variables['h0']
h = dat.variables['h']
h_min = dat.variables['h_min']
time = dat.variables['time']

#x_axis = np.linspace(-1,1,100)
#y_axis = np.linspace(-1,1,100)

fig, (ax0,ax1) = plt.subplots(1,2)

c = ax0.loglog(time[-1]-time,h_min)
c = ax0.loglog(time[-1]-time,time[-1]-time)
c = ax0.loglog(time[-1]-time,(time[-1]-time)**(1/5))
ax0.set_title('Logarithmic plot of time and h_min')
ax0.set_xlabel('time')
ax0.set_ylabel('h_min')
ax0.legend(['h_min','1','1/5'])
c = ax1.plot(x,h0,'r')
c = ax1.plot(x,h,'b')
ax1.set_title('shape of film')
ax1.set_xlabel('x')
ax1.set_ylabel('h')
plt.savefig('project3')
plt.show()
