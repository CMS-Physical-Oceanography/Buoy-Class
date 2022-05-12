import numpy as np
import pandas as pd 
from EddyClass import EddyVisc
from Coriolis import coriolis
class Simple:

    def __init__(self,h,PG,WS):
        # discritization
        self.dz = []
        self.Nt = 10000 # number of time-steps
        # Eddy Viscosity
        self.visc = EddyVisc() # Viscosity class
        self.Av = [] # array sent to Fortran
        # Data
        self.WndStrs = WS
        self.PresGrad = PG
        self.depth = h
        # Make sure data arrays are the same size
        self.check = self.datacheck()
    def datacheck(self):
        if len(self.WndStrs.i) == len(self.PresGrad.i):
            return 0 
        else:
            print('len(PresGrad)!=len(WndStrs)')
            return 1
    def makegrid(self,h,Nz,logscl=True):
        """This function creates a 1D spatial discritization 
        with depth using a logarithmic spacing away from the 
        boundaries if logscl=True. The finite difference scheme 
        doesnt always preform well with fine spatial discritization. 
        NOTE If a higher resolution is desired, the time-step, dt, 
        will need to be small and the model will take longer to reach a 
        steady state.
        INPUTS:
            -Nz = Number of vertical grid-points.
            -(logscl=True) = If true a logrithmic discritization 
                            will be assigned. If False dz will be 
                            depth uniform.
        REASSIGNS:
            -self.dz = Vertical discritization.
            """
        dz = h / (Nz-1)
        dz_ = np.zeros(Nz)
        idx= (Nz//2) + 1
        dz_[:idx] = -np.log(np.linspace(1.3,4.4,idx))*dz
        dz_[idx:] = -np.log(np.linspace(1.3,4.4,idx))[::-1][1:]*dz
        return dz_
            
    def makeAv(self,Nz,form):
        """
        This function calls
        Maybe define a lambda function for
        each form of Av that calls the 
        correct function from class: EddyVisc.

        THIS NEEDS TO BE CHANGED TO ASSIGN A FORM
        OF AV BASED ON THE INPUT: form.
        self.WndStrs.mag() for shear velocity dependent 
        forms of Av.
        -TH"""
        # check if PG and WS are the same length
        if self.check == 0:
            eddydict = {'bilincut':0,'const':1}
            self.dz = self.makegrid(self.depth,Nz)
            idx = eddydict[form]
            self.visc.dz = self.dz
            lendata = len(self.WndStrs.i)
            Av = np.zeros((Nz,lendata))
            for i in range(lendata):
                Av[:,i] = self.visc.Constant(.005)
            self.Av = Av # assign Av
        else:
            print('len(PresGrad)!=len(WndStrs)')
            return False 
    def tuneup(self,params):
        """
        This function will allow the user to tune the 
        model parameters using a python translation of
        the Fortran code in simple.f90.
        CONCEPT:
            - input wind stress and pressure gradient terms
              in a relativly small array and use them to
              force the python finite difference solver 
              which will be contained in Tuneup.py. This will
              use the Python code over Fortran to reduce the 
              setup speed assuming the input arrays are small.
              This will reassign the model parameters to the 
              inputted params with each call.
        -TH"""
        pass
        
    def pass2fortran(self):
        """
        This funciton passes class: Simple 
        attributes to Fortran. Data is read using 
        the setup.f90 module which uses fucntions
        from arryin.f90.
        -TH"""

        from ModelHelp import py2fort
        # Pass verticle discritization
        py2fort.dz2fortran(self.dz)
        # Pass time steps 
        py2fort.dt2fortran(self.Av,self.depth)
        # Pass Wind Stress log scale 
        py2fort.scale2fortran(self.Nt)
        # pass Av 
        py2fort.Av2fortran(self.Av)
        # pass wind stress data
        py2fort.data2fortran(self.WndStrs.i,self.WndStrs.j,
                                'csws.dat','asws.dat')
        # pass depth-average flow data
        py2fort.data2fortran(self.PresGrad.i,self.PresGrad.j,
                                'uda.dat','vda.dat')
    def results(self):
        """
        This function will read the model results from 
        Fortran binary and do something with them.
        -TH"""
        pass

    def CRUNCH(self):
        """This function compiles and executes the 
        Fortran90 files for the finite difference,
        and Stokes Drift profiles(soon as of 5/11/22).
       
        All discritization and data values must be assigned
        prior to calling CRUNCH().
        -TH"""
        import os
        import time 

        # Save Params and Data
        print('Passing run info to Fortran...')
        self.pass2fortran()
        # compile Fortran code
        print('Compiling...')
        os.system('cmd/c "gfortran -o test math.f90 Ocean.f90 arrayin.f90 savearry.f90 setup.f90 simple.f90"')
        start = time.time() # Start model run time
        print('Running Model...')
        os.system('cmd/c "test"') # Execute 
        end = time.time() # End model run time 

        runtime = (end-start)/60
        print('RUNTIME:',runtime,'minutes')

        self.results() # currently does nothing 
        # return runtime info
        return [self.Nt,len(self.dz),len(self.WndStrs.i),runtime]