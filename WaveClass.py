import numpy as np
from VectorClass import Vector2d

class Waves(Vector2d):
    def __init__(self, i,j,swh_=[], k_=[],h_ = 'Depth hasnt been set',cg_= [],T_ = []):
        super().__init__(i,j)
        self.traj = self.j
        self.swh = swh_
        self.k = k_
        self.T = T_
        self.cg = cg_
        self.h = h_
    
        

    def energyflux(self,swdens=1023,return_vector=True,zero_axis = 'y'):


        
        A = (swdens/16) * 9.81 
        Eflux = lambda swh,c,k : (  (A * swh**2 * c) * 
                                ( 1 + ((2* k* self.h) / ( np.sinh(2*k*self.h) )) )  )
        
        rad = np.deg2rad(self.j)
        E = np.array(list(map(Eflux,self.swh,self.i,self.k)))
    
        i = np.array(list(map(lambda r, theta : r*np.cos(theta),E,rad)))
        j = np.array(list(map(lambda r, theta : r*np.sin(theta),E,rad)))
        
        
        if return_vector == True:
            
            if zero_axis=='y':
                
                return Vector2d(j,i)
            
            elif zero_axis == 'x':
            
                return Vector2d(i,j)
            
            else:
                print('invalid input for zero_axis')
                pass
            
        elif return_vector == False:
            
            if zero_axis=='y':
                
                self.i = j
                self.j = i
                
            elif zero_axis == 'x':
                
                self.i = i
                self.j = j
                
            else:
                print('invalid input for zero_axis')
                pass
            
        else:

            print('invalid input to return_vector')
      

    def energy(self,swdens=1023):
        E = lambda swh :  ((swdens/8) * 9.81 * swh**2 )  
        out =   np.array(list(map(E,self.swh)))
        return out


    def getk(self):

        """
        ============================================================
        Uses the Newton-Raphson method to calculate the wavenumber
        k with an inital guess k0 being for shallow water waves where
        k0 = omega/sqrt(gh) and omega is the angular frequency 
        omega = 2pi/T. The function f is the dispersion relation for 
        surface gravity waves 0 = -omega + sqrt(gk*tanh(kh))
        This method is outlined in: 
        "A close approximation of wave dispersion relation for direct 
        calculation of wavelength in any coastal water depth" -Zai-Jin You

        This function does not return anything, instead it just assigns
        the calculated value to self.k
        ============================================================
                                                                 """
        
        # create omega and the initial guess of k defined in the doc string
        omega = (np.pi*2)/np.array(self.T)
        
        k = omega/np.sqrt(9.81*self.h) # make sure depth is set if error is thrown here 
    
        # create initial values using k0
        f = (9.81 * k * np.tanh(k*self.h)) - omega**2 
        dfdk = 9.81 * self.h * k * ( 1/(np.cosh(k*self.h)**2) ) + (9.81 * np.tanh(k*self.h))
        ii = np.any(abs(f) > 1e-10)
        
        while ii: # while ii == True
                                                                  
            k = k - (f/dfdk)
            f = (9.81 * k * np.tanh(k*self.h)) - omega**2
            dfdk = 9.81 * self.h * k * ( 1/(np.cosh(k*self.h)**2) ) + (9.81 * np.tanh(k*self.h))
            ii = np.any(abs(f) > 1e-10)
            
        self.k = k

    def getc(self):
        omega = (2*np.pi)/np.array(self.T)
        self.i = omega/self.k

    def getcg(self):

        Cg_ = lambda c,k,h :(  c*.5 * (1 + ((2*k*h)/np.sinh(2*k*h)))  )

        self.cg = Cg_(self.i,self.k,self.h)

        

    def stokesdrift(self,start = 'none',stop = 'none',res = 'none'):

        # The following three conditionals set the defult range of
        # stokes drift profiles to be calculated if they are not
        # provided by the user
        if stop == 'none':
            stop= len(self.swh)
        if start == 'none':
            start = 0
        if res  == 'none':
            res = int(self.h)
            
        # Makes sure that the depth is set
        if self.h ==  'Depth hasnt been set':
            print('set self.h')
            pass

        else:
            # initialize return and depth arrays
            xshelf = np.zeros((res,stop-start))
            yshelf = np.zeros((res,stop-start))
            z = np.linspace(-0.05,-self.h,res)

            # data slices used in calculations
            c = self.i[start:stop] # calculated with self.getc()
            theta = np.deg2rad(self.j[start:stop]) # 2/7/22: theta is degrees clockwise of +y
            swh = self.swh[start:stop]
            k_ = self.k[start:stop] # calculated with self.getk()

            # function to map z values through
            Sdrift_ = lambda z : (  np.cosh( 2*k*(z+self.h) )  ) # Only term in solution containing z
            
            for i in range(stop-start): # for each data point in [start:stop]

                k = k_[i] # needs to be decared for the Sdrift_ lambda function
                A =   ( (( 9.81*k*swh[i]**2) / ( 8*c[i] )) )
                # Cross-Shelf Stokes Drift
                xshelf[:,i] =( (A/ np.sinh( 2*k*self.h)) * np.cosh( 2*k*(z+self.h) ))*np.sin(theta[i]) 
                # Along-Shelf Stokes Drift
                yshelf[:,i] =( (A/ np.sinh( 2*k*self.h)) * np.cosh( 2*k*(z+self.h) )*np.cos(theta[i]) )
                
        return Vector2d(xshelf,yshelf) # returns a vector of stokes drift components

    def transport(self,zero_axis='y'):

        
        theta = np.deg2rad(self.j)
        Qw = lambda swh,cp : ( (9.81*swh**2)/
                                   (16*cp) )
        
        Q = np.array(Qw(self.swh,self.i))

        i = Q * np.sin(theta)
        j = Q * np.cos(theta)

        if zero_axis =='y':
            
            return Vector2d(j,i)
        
        elif zero_axis == 'x':

            return Vector2d([i,j])

        else:

            print('Invalid input to zero_axis')

   
    def waveforcing(self):
        
        Hsig = self.swh
        traj = np.deg2rad(self.j)
        
        return(Vector2d(Hsig**2 * np.sin(traj),Hsig**2 * np.cos(traj)))