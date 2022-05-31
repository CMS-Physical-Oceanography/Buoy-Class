import pandas as pd 
import numpy as np
"""
This file contains code used in ModelClass.py

Contains:
    -class: Py2Fort"""  

class py2fort:
    """
    Contains:
        -dz2fortran()
        -dt2fortran()
        -Av2fortran()
        -scale2fortran()
        -data2fortran()"""

    
    def dz2fortran(dz_):

        DZ = pd.DataFrame(dz_)
        DZ.to_csv('dz.dat',header=False,index=False)

    def dt2fortran(Av,h):  
        """
        This function calculates and saves the 
        timestep dt for each data point held in 
        class: Simple.
        INPUTS:
            -Av = Nz x lendata array of verticle 
                  eddy viscosity profiles.
            -h = Depth of the 1D model.
        SAVES:
            -dt = Time step for each viscosity profile."""
        sigma = .2 # arbitrary scale
        Nz = len(Av) # N vertical grid-points
        lendata = len(Av[0])
        dz = h / (Nz-1)

        # assign dt_ based on values of Av.
        # dt_ is a 1 x lendata array
        dt = np.array(list(map(lambda i : (sigma*dz**2)/max(Av[:,i]) 
                if max(Av[:,i]) > 0.002 else 50,range(lendata))))

        # assign all values of dt to zero where Av = 0
        # avoids errors in Fortran
        dt[list(map(lambda i : max(Av[:,i]==0),range(len(Av[0]))))] = 0

        # Save dt
        DT = pd.DataFrame(dt)
        DT.to_csv('dt.dat',header=False,index=False)

    def scale2fortran(Nt):
        wndscale_ = np.log(np.linspace(1.1,Nt,Nt)) 
        wndscale = np.zeros(len(wndscale_)+1)
        wndscale[0] = Nt
        wndscale[1:] = wndscale_/wndscale_[-1]
        WSC = pd.DataFrame(wndscale)
        WSC.to_csv('scale.dat',header=False,index=False)

    def Av2fortran(Av_):

        """
        This funcition outputs a 
        Nz x len(data) array of Eddy
        Viscosity profiles in a format
        that can be read by the setup.f90
        module using a function from arryin.f90."""

        Av = np.zeros(len(Av_[0])*len(Av_)+1)
        Av[0] = len(Av_) # Nz

        Av[1:] = Av_.flatten(order='F')    
        Av[np.isnan(Av)==True] = 0

        VISC = pd.DataFrame(Av)
        VISC.to_csv('Av.dat',header=False,index=False)

    def data2fortran(u,v,ufil,vfil):

        """
        This funciton creates a pressure 
        gradient array assuming the depth 
        averaged flow is in geostrophic 
        balence. Output is a len(data)
        complex array of (-vda + iuda)

        OUTPUT IS THE NEGATIVE PRESSURE GRAD
        WITHOUT THE CORIOLIS PARAMETER
        """

        v[np.isnan(v)==True] = 0
        u[np.isnan(u)==True] = 0

        v_ = np.zeros(len(v)+1)
        u_ = np.zeros(len(u)+1)

        v_[0] = len(v)
        u_[0] = len(u)

        u_[1:] = u
        v_[1:] = v

        U = pd.DataFrame(u_)
        V = pd.DataFrame(v_)

        U.to_csv(ufil,header=False,index=False)
        V.to_csv(vfil,header=False,index=False)

    def setup(simple):
        """
        This function inputs a Simple class object and passes
        its attruibutes to Fortran by saving them as .dat
        files."""

        print('Passing run info to Fortran...')
        # Pass verticle discritization
        py2fort.dz2fortran(simple.dz)
        # Pass time steps 
        py2fort.dt2fortran(simple.Av,simple.depth)
        # Pass Wind Stress log scale 
        py2fort.scale2fortran(simple.Nt)
        # pass Av 
        py2fort.Av2fortran(simple.Av)
        # pass wind stress data
        py2fort.data2fortran(simple.WS.i,simple.WS.j,
                                'csws.dat','asws.dat')
        # pass depth-average flow data
        py2fort.data2fortran(simple.PG.i,simple.PG.j,
                                'uda.dat','vda.dat')

class runfort:
    """
    This class contains functions used to run the
    Fortran model in Python."""

    
    def execute():

        import os
        from time import time

        # compile/run Fortran code and get run info
        print('Compiling...')
        os.system('cmd/c "gfortran -o test math.f90 Ocean.f90 arrayin.f90 savearry.f90 setup.f90 simple.f90"')
        start = time() # Start model run time
        print('Running Model...')
        os.system('cmd/c "test"') # Execute
        end = time() # End model run time 
        runtime = (end-start)/60
        print('RUNTIME:',runtime,'minutes')
        return runtime


    def clear():
        """
        This function deletes the data files passed to 
        Fortran in setup()."""
        
        from os import remove

        files = ['dz.dat','dt.dat','scale.dat','csws.dat','asws.dat','uda.dat','vda.dat']
        print('Clearing Inputs...')
        for fil in files:
            remove(fil)
        
