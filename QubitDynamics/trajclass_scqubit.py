'''
Trajectory Classes 
'''

class Traj:
    
    def __init__(self,system,r0,Pr,tvals,nps=1,tags=None,ori='plt',data=None,avN=1):
        self.system = system # Name of fSRK package containing system
        self.r0 = r0 # Array of initial conditions
        self.Pr = Pr # List of parameters
        self.tvals = tvals # Output time values 
        
        self.nps = nps # Number of integration steps per point in tvals
        self.ssd = len(r0) # State-space dimension        
        self.np = len(Pr) # Number of parameters
        
        self.tags = tags # Extra Information
        self.ori = ori # List of points: 'points' , List of trajectories for var. : 'plt'
        self.data = data
        
        self.avN = avN # AVerage of how many trajectories?
        
    def run(self):
        module = self.system
        rvals,ntvals = module.fsrk(r0=self.r0,t0=self.tvals[0],tf=self.tvals[-1],nout=len(self.tvals),nps=self.nps,p=self.Pr)
        self.data = rvals
        self.tvals = ntvals
        self.ori = 'plt'
        
    
class EnTraj:
    
    def __init__(self,Ntr,system,r0,Pr,tvals,nps=1,tags=None,ori='plt',data=None):
        self.Ntr = Ntr
        self.system = system # Name of fSRK package containing system
        self.r0 = r0 # Array of initial conditions
        self.Pr = Pr # List of parameters
        self.tvals = tvals # Output time values 
        
        self.nps = nps # Number of integration steps per point in tvals
        self.ssd = len(r0) # State-space dimension        
        self.np = len(Pr) # Number of parameters
        
        self.tags = tags # Extra Information
        self.ori = ori # List of points: 'points' , List of trajectories for var. : 'plt'
        self.data = data
        
    def run(self,timer=True):
        import xarray as xr
        from tqdm import tqdm
        from time import time
        
        if self.data == None:
            if timer == True:
                pbar = tqdm(total=self.Ntr-1)
                time1 = time()
            tmp_trj = Traj(self.system,self.r0,self.Pr,self.tvals,nps=self.nps)
            tmp_trj.run()
            tmp_foo =  xr.DataArray([tmp_trj.data],dims=['Ntr','ssd','Ntvl'],attrs={'r0':tmp_trj.r0,'Pr':tmp_trj.Pr,'tvals':[self.tvals[0],self.tvals[-1],len(self.tvals)]})
            tmp_data = tmp_foo
            for j in range(self.Ntr-1):
                if timer==True:
                    pbar.update(1)
                tmp_trj = Traj(self.system,self.r0,self.Pr,self.tvals,nps=self.nps)
                tmp_trj.run()
                tmp_foo =  xr.DataArray([tmp_trj.data],dims=['Ntr','ssd','Ntvl'],attrs={'r0':tmp_trj.r0,'Pr':tmp_trj.Pr,'tvals':[self.tvals[0],self.tvals[-1],len(self.tvals)]})
                tmp_ndata = xr.concat([tmp_data,tmp_foo],dim='Ntr')
                tmp_data = tmp_ndata
            self.data = tmp_data
            if timer==True:
                pbar.close()
                time2 = time()
                print('time taken:',time2-time1)
    
    def save(self,filename):
        self.data.to_netcdf(filename)
        
    def giveav(self):
        tedata = self.data
        return Traj(self.system,self.r0,self.Pr,self.tvals,data=tedata.sum(dim='Ntr')/self.Ntr,avN=self.Ntr)
    
    def givetraj(self,j):
        tedata = self.data 
        return Traj(self.system,self.r0,self.Pr,self.tvals,data=tedata[j])
    
    def give_avfunc(self,f):
        total = f(self.data[0])
        for i in range(1,len(self.data)):
            total = total + f(self.data[i])
        avfunc = [total/self.Ntr]
        return Traj(self.system,self.r0,self.Pr,self.tvals,data=avfunc, tags ='avfunc')

'''
Loading Functions
'''        
    
def load(filename,system):
    import xarray as xr
    import numpy as np
    dataarray = xr.open_dataarray(filename)
    x = dataarray.tvals
    loadeden = EnTraj(len(dataarray),system,dataarray.r0,dataarray.Pr,np.linspace(x[0],x[1],x[2]),data=dataarray)
    return loadeden

'''
Plotting Functions
'''

def pltcomp(traj,i=1,st=r'$z$'):
    import matplotlib.pyplot as plt
    data = traj.data
    fig = plt.figure(figsize =(5,4))
    plt.plot(traj.tvals,data[i])
    plt.xlabel(r'$t$',size=16)
    plt.ylabel(st,size=16)
    plt.grid(True)
    plt.show()
            
def pltphase(traj,i1=0,i2=1,st1=r'$y$',st2=r'$z$'):
    import matplotlib.pyplot as plt
    data = traj.data
    fig = plt.figure(figsize =(5,5))
    plt.plot(data[i1],data[i2])
    plt.xlabel(st1,size=16)
    plt.ylabel(st2,size=16)
    plt.grid(True)
    plt.show()
    
def pltfunc(traj,f,st):
    import matplotlib.pyplot as plt
    data = traj.data
    fig = plt.figure(figsize =(5,5))
    plt.plot(traj.tvals,f(data))
    plt.xlabel(r'$t$',size=16)
    plt.ylabel(st,size=16)
    plt.grid(True)
    plt.show()    
                
    
'''
Functions - Qubit 
'''

def purity(r):
    return (1 + r[0]**2 + r[1]**2 + r[2]**2)/2

def purity_sqb(r):
    return (1 + r[0]**2 + r[1]**2)/2
            
            