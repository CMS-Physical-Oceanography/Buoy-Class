import numpy as np
import pandas as pd
from .EddyClass import EddyVisc
import os

class Simple(EddyVisc):

    def __init__(self,depth=None,lat=None,WS=None,PG=None,buoy=None):
        super().__init__()
        self.WS,self.PG,self.depth,self.lat = self.initialize(depth,lat,WS,PG,buoy)
        self.check = self.datacheck()
        self.Nt = None

    def initialize(self,depth,lat,WS,PG,buoy):
        if buoy == None:
            return depth,WS,PG,lat

        else:
            wndstress = buoy.wind.windstress()
            flowda = buoy.currents.depth_average()
            return wndstress,flowda,buoy.depth,buoy.lat

    def datacheck(self):
        if len(self.WS.i) == len(self.PG.i):
            return True
        else:
            print('len(PG)!=len(WndStrs)')
            return False

    def makegrid(self,h,Nz,logscl=True):
        """
        This function creates a 1D spatial discritization 
        with depth using a logarithmic spacing away from the 
        boundaries if logscl=True. The finite difference scheme 
        doesnt always preform well with fine spatial discritization. 
        
        NOTE If a higher resolution is desired, the time-step, dt, 
        will need to be small and the model will take longer to reach a 
        steady state.
        INPUTS:
            -Nz = Number of vertical grid-points.
            -(logscl=True) = If true a logrithmic discritization 
            -               will be assigned. If False dz will be 
            -               depth uniform.
        RETURNS:
            -self.dz = Vertical discritization.
            """
        dz = h / (Nz-1)
        dz_ = np.zeros(Nz)
        idx= (Nz//2) + 1
        # dz_ = np.log(np.linspace(2.5,3.2,int(Nz)))*dz

        dz_[:idx] = np.log(np.linspace(1.3,4.4,idx))*dz
        dz_[idx:] = np.log(np.linspace(1.3,4.4,idx))[::-1][1:]*dz

        # dz_[:idx] = np.log(np.linspace(2.5,3.2,idx))*dz
        # dz_[idx:] = np.log(np.linspace(2.5,3.2,idx))[::-1][1:]*dz
        
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
        forms of Av."""

        # check if PG and WS are the same length
        if self.check == True:
            eddydict = {'bilincut':0,'const':1}
            self.dz = self.makegrid(self.depth,Nz)
            idx = eddydict[form]
            lendata = len(self.WS.i)
            self.Av = np.zeros((Nz,lendata))
            for i in range(lendata):
                self.Av[:,i] = self.Constant(.005)
        else:
            print('len(PG)!=len(WS)')
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
              inputted params with each call."""
        pass
        
    def pass2fortran(self):
        """
        This funciton passes class: Simple 
        attributes to Fortran. Data is read using 
        the setup.f90 module which uses fucntions
        from arryin.f90."""

        from .ModelHelp import py2fort
        py2fort.setup(self)
                 

    def clear_inputs(self,clean=True):
        from .ModelHelp import runfort
        if clean == True:
            runfort.clear()
        else:
            pass

    def results(self):
        """
        This function will read the model results from 
        Fortran binary and do something with them."""
        # from .CurrentClass import Currents 

        def replace_nans(arr):
            """
            when passed to Fortran nans are 
            replaced with 0. This replaces nans."""

            for col in range(len(arr[0])):
                if np.all(arr[:,col]==0):
                    arr[:,col] = np.nan
            return arr

        csflow = replace_nans(np.loadtxt('CSflow.bin'))
        asflow = replace_nans(np.loadtxt('ASflow.bin'))
        
        return csflow,asflow

    def run_model(self):
        """
        This function runs the Fortran code by calling 
        runfort.execute() and returns a run-info tuple."""

        from .ModelHelp import runfort
        runtime = runfort.execute() # compile/run Fortran code

        # return (Nt,Nz,len(data),runtime)
        return (self.Nt,len(self.dz),len(self.WS.i),runtime)

    def CRUNCH(self,clean=True):
        """
        This function compiles and executes the 
        Fortran90 files for the finite difference,
        and Stokes Drift profiles(soon as of 5/11/22).
       
        All discritization and data values must be assigned
        prior to calling CRUNCH()."""

        trace = os.getcwd()
        os.chdir('models')

        self.pass2fortran() # pass info to Fortran
        runinfo = self.run_model()
        self.clear_inputs(clean) # clear inputted data and params
        print('Done.')
        results = self.results()
        os.chdir(trace)
        return results,runinfo