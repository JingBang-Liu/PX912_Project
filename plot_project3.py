import netCDF4 as NC
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


# set text style
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text',usetext=True)

dat = NC.Dataset("project3.nc","r",format="NETCDF4")

#print(dat.__dict__)
#print(dat.variables['rho_charge_density'])
x = dat.variables['x']
h0 = dat.variables['h0']
h = dat.variables['h']
h_min = dat.variables['h_min']
time = dat.variables['time']

h_his = dat.variables['h_his']

h_min_size = h_min.size
T_temp = h_his.shape
T = T_temp[1]
print(T)
print('it took ',h_min_size,' iterations to reach ',h_min[h_min_size-2])


#x_axis = np.linspace(-1,1,100)
#y_axis = np.linspace(-1,1,100)

fig, (ax0,ax1) = plt.subplots(1,2)

c = ax0.loglog(time[-1]-time,h_min,'.',markersize=1,label=r'$h_{min}$')
c = ax0.loglog(time[-1]-time,time[-1]-time,'.',markersize=1,label='1')
c = ax0.loglog(time[-1]-time,(time[-1]-time)**(1/5),'.',markersize=1,label='$\\frac{1}{5}$')
c = ax0.loglog(time[-1]-time,(time[-1]-time)**(1/1000),'.',markersize=1,label='$\\frac{1}{100}$')
ax0.set_title('(a) Logarithmic plot of time and h_min')
ax0.set_xlabel('log(t_r-t)')
ax0.set_ylabel('log(h_min)')
#ax0.legend(['h_min','1','1/5'])
ax0.legend()
#ax0.text(0.5,0.05,'(a)',ha='center')
c = ax1.plot(x,h0,'r')
c = ax1.plot(x,h_his[:,1500],'g')
c = ax1.plot(x,h_his[:,1600],'b')
c = ax1.plot(x,h_his[:,1756],'c')
#c = ax1.plot(x,h,'b')
ax1.set_title('(b) shape of film')
ax1.set_xlabel('x')
ax1.set_ylabel('h')
#ax1.text(1.5,0.05,'(b)',ha='center')
plt.savefig('project3')
plt.show()
