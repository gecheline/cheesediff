import numpy as np 


class Cheese(object):

    def __init__(self, length, width, time = 0., c_salt = [], **kwargs):

        '''
        Initializes a cheese object with size, timestamp and measurements.

        Arguments:
        - length - length of cheese piece (in cm)
        - width - width of cheese piece (in cm)
        - time - the time measurement was done (in days since start of salting)
        - c_salt - the array of concentrations measured for the smaller 
        pieces cut up from the cheese piece
        
        '''

        self.length = length
        self.width = width
        self.time = time
        self.c = c_salt

        self.compute_centers_and_boundaries(**kwargs)


    def compute_centers_and_boundaries(self, nx=4, ny=3):

        '''
        Computes the centers and boundaries of the cut-up pieces

        Arguments:
        - nx - number of pieces in x direction
        - ny - number of pieces in y direction
        '''

        # check that the concentration array shape matches the provided number of pieces
        if len(self.c[0]) != nx*ny:
            raise ValueError('Concentration data doesn\'t match number of pieces: %i != %i x %i' 
                                % (len(self.c[0],nx,ny)))
        

        dl = self.length/nx
        dw = self.width/ny

        self._x_centers = np.array([dl*(i+0.5) for i in range(nx)])
        self._y_centers = np.array([dw*(j+0.5) for j in range(ny)])
        self._x_boundaries = np.array([i*dl for i in range(nx+1)])
        self._y_boundaries = np.array([j*dw for j in range(ny+1)])


    def add_brine_params(self, dl, dw, c_init = 1., c_brine = 1.):
        '''
        Extends the cheese centers, boundaries and measurements with brine data.

        Arguments:
        - dl - difference length of brine
        - dw - difference width of brine
        - c_init - initial salt concentration in brine
        - c_brine - measured concentration of brine in time=self.time
        '''
        # move center coordinates
        self.brine_dlength = dl
        self.brine_dwidth = dw 
        self.brine_c_salt = c_brine
        
        if not hasattr(self, 'brine_c0'):
            self.brine_c0 = c_init


class Diffusion(object):

    def __init__(self, length, width, nx, ny,  sigma):

        self.length = length
        self.width = width 
        self.xs = np.linspace(0,length,nx)
        self.ys = np.linspace(0,width,ny)
        self.dx = length/nx
        self.dy = width/ny
        self.sigma = 0.5


    def populate_init_c_box(self, dl_brine, dw_brine, c0_brine):
        # based on the values in the cheese and brine boundaries, populate init c
        self.dl_brine = dl_brine 
        self.dw_brine = dw_brine

        c0 = np.zeros((len(self.xs),len(self.ys)))
        argx = [np.argwhere(self.xs <= dl_brine).max(), np.argwhere(self.xs >= self.length - dl_brine).min()]
        argy = [np.argwhere(self.ys <= dw_brine).max(), np.argwhere(self.ys >= self.width - dw_brine).min()]

        c0[:argx[0],:] = c0_brine
        c0[argx[1]:,:] = c0_brine
        c0[:,:argy[0]] = c0_brine
        c0[:,argy[1]:] = c0_brine

        mask = np.zeros(c0.shape)
        mask[c0==c0_brine] = 1
        self.c0 = c0
        self.brine_mask = mask.astype(bool)

    def diffuse(self, D, nt=1):
        
        self.diffusion_constant = D
        self.dt = self.sigma * (self.dx**2*self.dy**2) /(2*D*(self.dx**2+self.dy**2))
        
        c_init = self.c0
        c_new = c_init.copy()

        for i in range(nt):
            c_new[1:-1,1:-1] = c_init[1:-1,1:-1] + D*self.dt/self.dx**2*(c_init[2:,1:-1]-2*c_init[1:-1,1:-1]+c_init[0:-2,1:-1]) + D*self.dt/self.dy**2*(c_init[1:-1,2:]-2*c_init[1:-1,1:-1]+c_init[1:-1,0:-2]) 
            c_new[0,:] = c_new[1,:]
            c_new[-1,:] = c_new[-2,:]
            c_new[:,0] = c_new[:,1]
            c_new[:,-1] = c_new[:,-2]

            c_init = c_new.copy()
        
        self.c_diffusion = c_new


    def compute_averages(self, cheese_boundaries, map=[[9,10,11,12],[5,6,7,8],[1,2,3,4]]):
        # taking the boundaries compute averages for cheese slices and brine and return as array
        if not hasattr(self, 'c_diffusion'):
            raise ValueError('No computed diffusion concentrations! Run .diffuse() first!')
        
        #let's separate cheese and brine
        cheese_xs = self.xs[~self.brine_mask[int(len(self.brine_mask)/2)]]
        cheese_ys = self.ys[~self.brine_mask.T[int(len(self.brine_mask.T)/2)]]

        brine = self.c_diffusion[self.brine_mask]
        cheese = self.c_diffusion[~self.brine_mask].reshape((len(cheese_xs),len(cheese_ys)))

        c_brine = np.average(brine)

        # figure out the format of the cheese boundaries and how to compute the average
        c_cheese = []

        return None

    def animate(self, nt=10):
        return None

    def plot_c_grid(self):
        return None

    def plot_c_averages(self):
        return None

class Fitter(object):

    def __init__(self, cheese_params={}, brine_params={}, diffusion_params={}):

        # take parameters for the cheese and brine data and the simulation
        return None

    def fit_model(self):
        return None


    def minimize(self):
        return None

    
    def sample(self):
        return None


    def plot_results(self):
        return None
