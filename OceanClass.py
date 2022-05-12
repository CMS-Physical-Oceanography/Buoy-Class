"""
OceanClass contains classes for dealing with wind, wave, and ADCP current measurements.
Since wind and wave directional data is typically reported in degrees clockwise from the
positive y axis (North) functions assume that sin(angle) is the x component (East) and 
cos(angle) is the y component

Data pulled from NOAA's NDBC database will need to be inverted from direction of origin to
direction of movement using the self.invert() function in the Vector2d parent class"""
import numpy as np
from VectorClass import Vector2d
from WindClass import Wind
from WaveClass import Waves
from CurrentClass import Currents
from EddyClass import EddyVisc
from ModelClass import Simple
from BuoyClass import Buoy
from Coriolis import coriolis


# class Waves(Vector2d):
#     def __init__(self, data, swh_=[], k_=[],h_ = 'Depth hasnt been set',cg_= [],T_ = []):
#         super().__init__(data)
#         self.traj = self.j
#         self.swh = swh_
#         self.k = k_
#         self.T = T_
#         self.cg = cg_
#         self.h = h_
    
        

#     def energyflux(self,swdens=1023,return_vector=True,zero_axis = 'y'):


        
#         A = (swdens/16) * 9.81 
#         Eflux = lambda swh,c,k : (  (A * swh**2 * c) * 
#                                 ( 1 + ((2* k* self.h) / ( np.sinh(2*k*self.h) )) )  )
        
#         rad = np.deg2rad(self.j)
#         E = np.array(list(map(Eflux,self.swh,self.i,self.k)))
    
#         i = np.array(list(map(lambda r, theta : r*np.cos(theta),E,rad)))
#         j = np.array(list(map(lambda r, theta : r*np.sin(theta),E,rad)))
        
        
#         if return_vector == True:
            
#             if zero_axis=='y':
                
#                 return Vector2d([j,i])
            
#             elif zero_axis == 'x':
            
#                 return Vector2d([i,j])
            
#             else:
#                 print('invalid input for zero_axis')
#                 pass
            
#         elif return_vector == False:
            
#             if zero_axis=='y':
                
#                 self.i = j
#                 self.j = i
                
#             elif zero_axis == 'x':
                
#                 self.i = i
#                 self.j = j
                
#             else:
#                 print('invalid input for zero_axis')
#                 pass
            
#         else:

#             print('invalid input to return_vector')
      

#     def energy(self,swdens=1023):
#         E = lambda swh :  ((swdens/8) * 9.81 * swh**2 )  
#         out =   np.array(list(map(E,self.swh)))
#         return out


#     def getk(self):

#         """
#         ============================================================
#         Uses the Newton-Raphson method to calculate the wavenumber
#         k with an inital guess k0 being for shallow water waves where
#         k0 = omega/sqrt(gh) and omega is the angular frequency 
#         omega = 2pi/T. The function f is the dispersion relation for 
#         surface gravity waves 0 = -omega + sqrt(gk*tanh(kh))
#         This method is outlined in: 
#         "A close approximation of wave dispersion relation for direct 
#         calculation of wavelength in any coastal water depth" -Zai-Jin You

#         This function does not return anything, instead it just assigns
#         the calculated value to self.k
#         ============================================================
#                                                                  """
        
#         # create omega and the initial guess of k defined in the doc string
#         omega = (np.pi*2)/np.array(self.T)
        
#         k = omega/np.sqrt(9.81*self.h) # make sure depth is set if error is thrown here 
    
#         # create initial values using k0
#         f = (9.81 * k * np.tanh(k*self.h)) - omega**2 
#         dfdk = 9.81 * self.h * k * ( 1/(np.cosh(k*self.h)**2) ) + (9.81 * np.tanh(k*self.h))
#         ii = np.any(abs(f) > 1e-10)
        
#         while ii: # while ii == True
                                                                  
#             k = k - (f/dfdk)
#             f = (9.81 * k * np.tanh(k*self.h)) - omega**2
#             dfdk = 9.81 * self.h * k * ( 1/(np.cosh(k*self.h)**2) ) + (9.81 * np.tanh(k*self.h))
#             ii = np.any(abs(f) > 1e-10)
            
#         self.k = k

#     def getc(self):
#         omega = (2*np.pi)/np.array(self.T)
#         self.i = omega/self.k

#     def getcg(self):

#         Cg_ = lambda c,k,h :(  c*.5 * (1 + ((2*k*h)/np.sinh(2*k*h)))  )

#         self.cg = Cg_(self.i,self.k,self.h)

        

#     def stokesdrift(self,start = 'none',stop = 'none',res = 'none'):

#         # The following three conditionals set the defult range of
#         # stokes drift profiles to be calculated if they are not
#         # provided by the user
#         if stop == 'none':
#             stop= len(self.swh)
#         if start == 'none':
#             start = 0
#         if res  == 'none':
#             res = int(self.h)
            
#         # Makes sure that the depth is set
#         if self.h ==  'Depth hasnt been set':
#             print('set self.h')
#             pass

#         else:
#             # initialize return and depth arrays
#             xshelf = np.zeros((res,stop-start))
#             yshelf = np.zeros((res,stop-start))
#             z = np.linspace(-0.05,-self.h,res)

#             # data slices used in calculations
#             c = self.i[start:stop] # calculated with self.getc()
#             theta = np.deg2rad(self.j[start:stop]) # 2/7/22: theta is degrees clockwise of +y
#             swh = self.swh[start:stop]
#             k_ = self.k[start:stop] # calculated with self.getk()

#             # function to map z values through
#             Sdrift_ = lambda z : (  np.cosh( 2*k*(z+self.h) )  ) # Only term in solution containing z
            
#             for i in range(stop-start): # for each data point in [start:stop]

#                 k = k_[i] # needs to be decared for the Sdrift_ lambda function
#                 A =   ( (( 9.81*k*swh[i]**2) / ( 8*c[i] )) )
#                 # Cross-Shelf Stokes Drift
#                 xshelf[:,i] =( (A/ np.sinh( 2*k*self.h)) * np.cosh( 2*k*(z+self.h) ))*np.sin(theta[i]) 
#                 # Along-Shelf Stokes Drift
#                 yshelf[:,i] =( (A/ np.sinh( 2*k*self.h)) * np.cosh( 2*k*(z+self.h) )*np.cos(theta[i]) )
                
#         return Vector2d([xshelf,yshelf]) # returns a vector of stokes drift components

#     def transport(self,zero_axis='y'):

        
#         theta = np.deg2rad(self.j)
#         Qw = lambda swh,cp : ( (9.81*swh**2)/
#                                    (16*cp) )
        
#         Q = np.array(Qw(self.swh,self.i))

#         i = Q * np.sin(theta)
#         j = Q * np.cos(theta)

#         if zero_axis =='y':
            
#             return Vector2d([j,i])
        
#         elif zero_axis == 'x':

#             return Vector2d([i,j])

#         else:

#             print('Invalid input to zero_axis')

   
#     def waveforcing(self):
        
#         Hsig = self.swh
#         traj = np.deg2rad(self.j)
        
#         return(Vector2d([Hsig**2 * np.sin(traj),Hsig**2 * np.cos(traj)]))
                


# class Wind(Vector2d):

#     """
#     Wind  inputs wind speed as self.i and direction as 
#     self.j
 
#     Like its parent class Vector2d, built in functions in 
#     Wind assume self.i is the radial or horizontal 
#     component and self.j is the angular or vertical component
#     of the wind vector. To work with typical meterological data
#     angles are taken to be clockwise from the positive y axis
    
    
                                                            
#     __init__():

#     self.i
#     self.j
#     self.cd # drag coeff.
    
#     1) windstress(self, drag = 1.4E-3, airdens = 1.22, 
#               return_vector = True,zero_axis = 'y',input_coordinates = 'polar)
              
#        self.windstress() assumes that that inputed windspeeds and drag coefficients 
#        are adjusted to the 10 meter wind speed equivalents. If data is collected at 
#        a different height call self.cdnlp before

#                                                                """

#     def __init__(self,data):
        
#         super().__init__(data)
#         self.cd = np.ones(len(self.i))*1.15e-3 


#     def cdnlp(self,heights):

#         """
#         =======================================================
#         MODIFIED FROM:
#         https://github.com/pyoceans/python-airsea/blob/master/airsea/windstress.py
        
#         The python-airsea toolbox was translated from the 
#         matlab version.  

#         This functions calculates U10 and drag coefficient cd 
#         at different heights and windspeeds for the values in 
#         self.i as outlined in Large & Pond 1981. 

#         inputs:
        
#         self

#         heights = vertical distance from the mesurement height
#                   to the surface roughness(Large & Pond) in meters.

#         returns: Nothing

#         re-assigns: 
        
#         self.i = u10
        
#         self.cd = cd

#         =======================================================
#         """
#         sp = self.i
        
#         tol = 0.001  # minimum convergence 
                  
           
#         z = np.array(heights)
#         a = np.log(z / 10.) / .4  # Log-layer correction factor.
#         u10o = np.zeros(sp.shape)
#         cd = 1.15e-3 * np.ones(sp.shape)
#         u10 =  sp / (1 + a * np.sqrt(cd))
#         ii = np.abs(u10 - u10o) > tol
        
#         while np.any(ii):
#             u10o = u10
#             cd = (4.9e-4 + 6.5e-5 * u10o)  # Compute cd(u10).
#             cd[u10o < 10.15385] = 1.15e-3
#             u10 = sp / (1 + a * np.sqrt(cd))  # Next iteration.
#             # Keep going until iteration converges.
#             ii = np.abs(u10 - u10o) >tol

#         self.i = u10
#         self.cd = cd
        

#     def Ekmantransport(self,swdens = 1023,lat = 33.5):

#         wstrs = self.windstress(return_vector= True)
        
#         cor = coriolis(lat)
#         Uekx = wstrs.j /(swdens*cor)
#         Ueky = wstrs.i / (swdens*cor)

#         return Vector2d([Uekx,Ueky])



#     def windstress(self, airdens = 1.22,real=True):
        
#         """

#         MAKE SURE TO CALL SELF.CDNLP OR DRAG
#         WILL BE SET TO A DEFULT VALUE OF 
#         1.15e-3.
#         """
 
#         drag = self.cd
       
#         rad = np.deg2rad(self.j)
#         wndstress = drag*airdens*self.i**2
#         j = wndstress*np.cos(rad)
#         i = wndstress*np.sin(rad)

#         if real==False:

#             out = np.zeros(len(self.i),dtype=np.complex_)
#             out.real = i
#             out.imag = j
#             return out

#         else:

#             return Vector2d([j,i])


    
# class Currents(Vector2d):
    
#     """
#     Currents Class expects input data to have dimentions depth by time
#     or depth by spatial dim. Angles are reported from the positive y axis.
#     """
    
#     def __init__(self,csflow,asflow,depths_ = np.nan):
        
#         self.i = csflow
#         self.j = asflow
#         self.z = depths_
        
        
#     def rot_angles(self,theta_prime):
        
        
        
#         u_0 = self.i
#         v_0 = self.j
        
        
#         # calculates radial comp
#         r_0 = np.sqrt(u_0**2 + v_0**2)
#         # calculates angle from the positive y
#         theta_0 = np.arctan2(u_0,v_0)
    
#         # maps angles to new axis
#         theta_f = theta_0 - np.deg2rad(theta_prime)
       
#         # calculates components for new coordinate system
#         self.i = r_0 * np.sin(theta_f) # cross-shelf
#         self.j = r_0 * np.cos(theta_f) # along-shelf
        


#     def transport(self):



#         tx = range(len(self.i[0]))
#         ty = range(len(self.j[0]))
        
#         cross_shelf = lambda t : np.nansum(self.i[:,t])
#         along_shelf = lambda t : np.nansum(self.j[:,t])

#         cs_transports = np.array(list(map(cross_shelf,tx)))
#         as_transports = np.array(list(map(along_shelf,ty)))

#         return Vector2d([cs_transports,as_transports])



#     def depth_average(self,real=True):



#         # tx = range(len(self.i[0]))
#         # ty = range(len(self.j[0]))
        
#         # cross_shelf = lambda t : np.nanmean(self.i[:,t])
#         # along_shelf = lambda t : np.nanmean(self.j[:,t])

#         # cs_da = np.array(list(map(cross_shelf,tx)))
#         # as_da = np.array(list(map(along_shelf,ty)))

#         cs_da = np.nanmean(self.i,axis=0)
#         as_da = np.nanmean(self.j,axis=0)


#         if real == False:

#             out = np.zeros(self.i.shape[1],dtype=np.complex_)
#             out.real = cs_da
#             out.imag = as_da
#             return out 

#         return Vector2d([cs_da,as_da]) 

# # class Ekman_Spiral:
# #     def __init__(self,A,delta):
# #         self.u = list(map(lambda z: A*np.exp(z/delta)*
# #                      np.sin((np.pi/4)-(z/delta)),range(-h,1)))
# #         self.v = list(map(lambda z: A*np.exp(z/delta)*
# #                      np.cos((np.pi/4)-(z/delta)),range(-h,1)))