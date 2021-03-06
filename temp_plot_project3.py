import netCDF4 as NC
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


# set text style
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text',usetext=True)

dat = NC.Dataset("project3_temp.nc","r",format="NETCDF4")

#print(dat.__dict__)
#print(dat.variables['rho_charge_density'])
x = dat.variables['x']
h0 = dat.variables['h0']
h = dat.variables['h']
h_min = dat.variables['h_min']
time = dat.variables['time']
h_his = dat.variables['h_his']

T_temp = h_his.shape
T = T_temp[1]
print(T)

#x_axis = np.linspace(-1,1,100)
#y_axis = np.linspace(-1,1,100)

fig, (ax0,ax1) = plt.subplots(1,2)

c = ax0.loglog(time[-1]-time,h_min,'.',label=r'$h_{min}$')
c = ax0.loglog(time[-1]-time,time[-1]-time,label='1')
c = ax0.loglog(time[-1]-time,(time[-1]-time)**(1/5),label='$\\frac{1}{5}$')
ax0.set_title('Logarithmic plot of time and h_min')
ax0.set_xlabel('log(t_r-t)')
ax0.set_ylabel('log(h_min)')
#ax0.legend(['h_min','1','1/5'])
ax0.legend()
c = ax1.plot(x,h0,'r')
#c = ax1.plot(x,h_his[:,1:10:70],'b')
c = ax1.plot(x,h_his[:,71],'b')
c = ax1.plot(x,h_his[:,72],'g')
c = ax1.plot(x,h_his[:,70],'c')
ax1.set_title('shape of film')
ax1.set_xlabel('x')
ax1.set_ylabel('h')
plt.savefig('project3')
plt.show()
