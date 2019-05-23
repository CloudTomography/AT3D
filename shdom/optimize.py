"""
Optimization and related objects to monitor and log the optimization proccess.
"""
import shdom
import numpy as np
import core, time, os
from scipy.optimize import minimize
from shdom import GridData
import dill as pickle
import cv2

    
class SummaryWriter(object):
    """
    A wrapper for tensorboardX summarywriter with some basic summary writing implementation.
    This wrapper enables logging of images, error measures and loss with pre-determined temporal intervals into tensorboard.
    To see the summary of this run (and comparisons to all subdirectories):
        tensorboard --logdir LOGDIR 
    """
    def __init__(self, log_dir):
        
        from tensorboardX import SummaryWriter
        
        self._tf_writer = SummaryWriter(log_dir)
        self._ground_truth_parameters = None
        self._callback_fns = []
        self._optimizer = None
        
        
    def attach_optimizer(self, optimizer):
        self._optimizer = optimizer
    
    def monitor_images(self, acquired_images, ckpt_period=-1):
        """
        Monitor the synthetic images and compare to the acquired images
        
        Parameters
        ----------
        acquired_images: list
            List of aquired images will be logged once onto tensorboard for comparison with the current state.
        ckpt_period: float
           time [seconds] between updates. setting ckpt_period=-1 will log at every iteration.
        """
        self._image_ckpt_period = ckpt_period
        self._image_ckpt_time = time.time()
        self._image_vmax = []
        for i, image in enumerate(acquired_images):
            self._image_vmax.append(image.max() * 1.25)
            self.tf_writer.add_image(tag='Acquiered Image{}'.format(i), 
                                     img_tensor=(image / self._image_vmax[i]),
                                     dataformats='HW')        
        self._callback_fns.append(self.images)
    
    
    def images(self):
        """Callback function the is called every optimizer iteration image monitoring is set."""
        time_passed = time.time() - self._image_ckpt_time 
        if time_passed > self._image_ckpt_period:
            self._image_ckpt_time  = time.time()
            for i, image in enumerate(self.optimizer.images):
                self.tf_writer.add_image(tag='Retrieval Image{}'.format(i), 
                                         img_tensor=(image / self._image_vmax[i]), 
                                         global_step=self.optimizer.iteration,
                                         dataformats='HW')
        
        
    def monitor_loss(self):
        """Monitor the loss at every iteration."""
        self._callback_fns.append(self.loss)
    
    
    def loss(self):
        """Callback function the is called every optimizer iteration loss monitoring is set."""
        self.tf_writer.add_scalar('loss', self.optimizer.loss, self.optimizer.iteration)
        
        
    def monitor_parameter_error(self, ground_truth_params, ckpt_period=-1):
        """
        Monitor relative and overall mass error (epsilon, delta) as defined by:
          Amit Aides et al, "Multi sky-view 3D aerosol distribution recovery".
        
        Parameters
        ----------
        ground_truth_params: list
            A list of the ground truth parameters to compare current state with.
        ckpt_period: float
           time [seconds] between updates. setting ckpt_period=-1 will log at every iteration.
        """
        self._param_ckpt_period = ckpt_period
        self._param_ckpt_time = time.time()        
        self._ground_truth_params = ground_truth_params
        self._callback_fns.append(self.param_error)
        
        
    def param_error(self):
        """Callback function the is called every optimizer iteration parameter error monitoring is set."""
        time_passed = time.time() - self._image_ckpt_time 
        if time_passed > self._param_ckpt_period:
            for est_param, gt_param in zip(self.optimizer.parameters, self._ground_truth_params):
                delta = (np.linalg.norm(est_param.data.ravel(),1) - np.linalg.norm(gt_param.data.ravel(),1)) / np.linalg.norm(gt_param.data.ravel(),1)
                epsilon = np.linalg.norm((est_param - gt_param).data.ravel(),1) / np.linalg.norm(gt_param.data.ravel(),1)
                self.tf_writer.add_scalar('delta', delta, self.optimizer.iteration)
                self.tf_writer.add_scalar('epsilon', epsilon, self.optimizer.iteration)           
    
    @property
    def callback_fns(self):
        return self._callback_fns
    
    @property
    def optimizer(self):
        return self._optimizer
    
    @property
    def tf_writer(self):
        return self._tf_writer
    
    
        

class SpaceCarver(object):
    """
    SpaceCarver object recovers the convex hull of the cloud based on multi-view sensor geometry.
    
    Parameters
    ----------
    measurements: shdom.Measurements
        A measurements object storing the images and sensor geometry
    """
    def __init__(self, measurements):
    
        self._rte_solver = shdom.RteSolver(shdom.SceneParameters(), shdom.NumericalParameters())
        
        self._measurements = measurements
        
        if measurements.sensors.type == 'SensorArray':
            self._sensors = measurements.sensors.sensor_list
        else:
            self._sensors = [measurements.sensors]
        self._images = measurements.images
        
        
    def carve(self, grid, agreement=0.75):
        """
        Carves out the cloud geometry on the grid. 
        A threshold on the radiances is found by a cv2 adaptive image threshold.
        
        Parameters
        ----------
        grid: shdom.Grid
            A grid object.
        agreement: the precentage of pixels that should agree on a cloudy voxels to set it to True in the mask
        
        Returns
        -------
        mask: shdom.GridData object
            A boolean mask with True marking cloudy voxels and False marking non-cloud region.
        """
           
        self._rte_solver.set_grid(grid)
        volume = np.zeros((grid.nx, grid.ny, grid.nz))
        
        for sensor, image in zip(self._sensors, self._images):
            uint8 = np.uint8(255*image/image.max())
            threshold = cv2.threshold(uint8, 0, 1, cv2.THRESH_BINARY + cv2.THRESH_OTSU)[1]           

            radiance = image.ravel(order='F')
            sensor = sensor[threshold.ravel(order='F') == 1]
            
            carved_volume = core.space_carve(
                nx=grid.nx,
                ny=grid.ny,
                nz=grid.nz,
                npts=self._rte_solver._npts,
                ncells=self._rte_solver._ncells,
                gridptr=self._rte_solver._gridptr,
                neighptr=self._rte_solver._neighptr,
                treeptr=self._rte_solver._treeptr,
                cellflags=self._rte_solver._cellflags,
                bcflag=self._rte_solver._bcflag,
                ipflag=self._rte_solver._ipflag,
                xgrid=self._rte_solver._xgrid,
                ygrid=self._rte_solver._ygrid,
                zgrid=self._rte_solver._zgrid,
                gridpos=self._rte_solver._gridpos,
                camx=sensor.x,
                camy=sensor.y,
                camz=sensor.z,
                cammu=sensor.mu,
                camphi=sensor.phi,
                npix=sensor.npix,
            )
            volume += carved_volume.reshape(grid.nx, grid.ny, grid.nz)
        
        volume = volume * 1.0 / len(self._images)
        mask = GridData(grid, volume > agreement) 
        return mask
    
    
    @property
    def grid(self):
        return self._grid


class Optimizer(object):
    """
    The Optimizer class takes care of the under the hood of the optimization process. 
    To run the optimization the following methods should be called:
       [required] optimizer.set_measurements()
       [required] optimizer.set_rte_solver()
       [optional] optimizer.set_cloud_mask()
       [required] optimizer.add_parameter(parameter) 
       [optional] optimizer.set_known_medium(air)
       [optional] optimizer.set_writer(writer) 
    
    Notes
    -----
    Currently only Extinction optimization is supported.
    """
    def __init__(self):
        self._parameters = []
        self._known_medium = None
        self._rte_solver = None
        self._measurements = None
        self._cloud_mask = None
        self._gradient_mask = None
        self._num_parameters = 0
        self._extinction_dependency = False
        self._albedo_dependency = False
        self._phase_dependency = False 
        self._writer = None
        self._images = None
        self._iteration = 0
        self._ckpt_time = time.time()
        self._ckpt_period = None
        self._loss = None
        
    def set_measurements(self, measurements):
        """
        Set the measurements (data-fit constraints)
        
        Parameters
        ----------
        measurements: shdom.Measurements
            A measurements object storing the images and sensor geometry
        """
        self._measurements = measurements
        
        
    def set_rte_solver(self, rte_solver):
        """
        Set the RteSolver for the SHDOM iterations.
        
        Parameters
        ----------
        rte_solver: shdom.RteSolver object
            The RteSolver initialized to the medium (RteSolver.init_medium() method).
        """
        self._rte_solver = rte_solver

    def set_cloud_mask(self, cloud_mask):
        """
        Set the cloud mask for the gradient descent.
        
        Parameters
        ----------
        cloud_mask: shdom.GridData object
            A boolean mask with True making cloudy voxels and False marking non-cloud region.
        
        Notes
        -----
        A None cloud mask means no mask
        """
        self._cloud_mask = cloud_mask  
    
    def set_known_medium(self, medium):
        """
        The known medium (e.g. known molecular scattering) is not optimized for but is accounted for in the forward model.
        
        Parameters
        ----------
        medium: shdom.Medium
            A medium object.
        """
        self._known_medium = medium


    def set_writer(self, writer):
        """
        Set a log writer to upload summaries into tensorboard.
        
        Parameters
        ----------
        writer: shdom.SummaryWriter
            Wrapper for the tensorboardX summary writer.
        """
        self._writer = writer
        self._writer.attach_optimizer(self)
        
        
    def add_parameter(self, parameter):
        """
        Add a Parameter to optimize for. 
        
        Parameters
        ----------
        parameter: shdom.Parameter
            A parameter for the optimization.
        
        Notes
        -----
        Currently only Extinction is supported.
        """
        if parameter.__class__ is not shdom.parameters.Extinction:
            raise NotImplementedError('Only Extinction optimization is supported.')
        
        self._parameters.append(parameter)
        self._extinction_dependency = self.extinction_dependency or parameter.extinction_dependency
        self._albedo_dependency = self.albedo_dependency or parameter.albedo_dependency
        self._phase_dependency = self.phase_dependency or parameter.phase_dependency
        
        
    def update_rte_solver(self, state):
        """
        An internal function to update the RteSolver to the current state.
        
        state: np.array(dtype=float64)
            The state is a numpy array of all the parameters concatenated. This is required by the scipy optimizer.
        """ 
        param_data = np.split(state, np.cumsum(self.num_parameters[:-1]))
        
        for i, param in enumerate(self.parameters):
            param.set_data(param_data[i])
            
            if param.extinction_dependency:
                extinction = param.get_extinction()
                if self.known_medium:
                    extinction += self.known_medium.extinction
                self.rte_solver.set_extinction(extinction)
            
            if param.albedo_dependency:
                raise NotImplementedError                
            
            if param.phase_dependency:
                raise NotImplementedError
            
            
    def extinction_gradient_cost(self):
        """
        Gradient and cost computation by ray tracing the domain (see shdomsub4.f).
        Additionally the synthetic images produces by the sensors at the current state are saved (for logging/visualization)
        
        Returns
        -------
        gradient: np.array(dtype=float64)
            The gradient of the objective function at the current state.
        cost: float
            The objective function (cost) evaluated at the current state
        
        Notes
        -----
        Currently the forward SHDOM solution is invoked at every gradient computation. 
        """ 
        
        self.rte_solver.solve(maxiter=100, verbose=False)   
        gradient, cost, images = core.ext_gradient(
            nx=self.rte_solver._nx,
            ny=self.rte_solver._ny,
            nz=self.rte_solver._nz,
            bcflag=self.rte_solver._bcflag,
            ipflag=self.rte_solver._ipflag,   
            npts=self.rte_solver._npts,
            nbpts=self._rte_solver._nbpts,
            ncells=self.rte_solver._ncells,
            nbcells=self.rte_solver._nbcells,
            ml=self.rte_solver._ml,
            mm=self.rte_solver._mm,
            ncs=self.rte_solver._ncs,
            nlm=self.rte_solver._nlm,
            numphase=self.rte_solver._pa.numphase,
            nmu=self.rte_solver._nmu,
            nphi0max=self.rte_solver._nphi0max,
            nphi0=self.rte_solver._nphi0,
            maxnbc=self.rte_solver._maxnbc,
            ntoppts=self.rte_solver._ntoppts,
            nbotpts=self.rte_solver._nbotpts,
            nsfcpar=self.rte_solver._nsfcpar,
            gridptr=self.rte_solver._gridptr,
            neighptr=self.rte_solver._neighptr,
            treeptr=self.rte_solver._treeptr,             
            shptr=self.rte_solver._shptr,
            bcptr=self.rte_solver._bcptr,
            cellflags=self.rte_solver._cellflags,
            iphase=self.rte_solver._iphase,
            deltam=self.rte_solver._deltam,
            solarmu=self.rte_solver._solarmu,
            solaraz=self.rte_solver._solaraz,
            gndtemp=self.rte_solver._gndtemp,
            gndalbedo=self.rte_solver._gndalbedo,
            skyrad=self.rte_solver._skyrad,
            waveno=self.rte_solver._waveno,
            wavelen=self.rte_solver._wavelen,
            mu=self.rte_solver._mu,
            phi=self.rte_solver._phi.reshape(self.rte_solver._nmu, -1),
            wtdo=self.rte_solver._wtdo.reshape(self.rte_solver._nmu, -1),
            xgrid=self.rte_solver._xgrid,
            ygrid=self.rte_solver._ygrid,
            zgrid=self.rte_solver._zgrid,
            gridpos=self.rte_solver._gridpos,
            sfcgridparms=self.rte_solver._sfcgridparms,
            bcrad=self.rte_solver._bcrad,
            extinct=self.rte_solver._extinct,
            albedo=self.rte_solver._albedo,
            legen=self.rte_solver._legen.reshape(self.rte_solver._nleg+1, -1),            
            dirflux=self.rte_solver._dirflux,
            fluxes=self.rte_solver._fluxes,
            source=self.rte_solver._source,
            camx=self.sensors.x,
            camy=self.sensors.y,
            camz=self.sensors.z,
            cammu=self.sensors.mu,
            camphi=self.sensors.phi,
            npix=np.sum(self.sensors.npix),          
            srctype=self.rte_solver._srctype,
            sfctype=self.rte_solver._sfctype,
            units=self.rte_solver._units,
            measurements=self.radiances,
            rshptr=self.rte_solver._rshptr,
            radiance=self.rte_solver._radiance
        )
        
        # Split images into different sensor images
        if self.sensors.type == 'SensorArray':
            image_list = np.split(images, np.cumsum(self.sensors.npix[:-1]))
            sensor_list = self.sensors.sensor_list
        else:
            image_list = [images]
            sensor_list = [self.sensors]
            
        for i in range(len(image_list)):
            if hasattr(sensor_list[i], 'resolution'):
                image_list[i] = image_list[i].reshape(sensor_list[i].resolution, order='F')
        self._images = image_list
        
        return gradient, cost
    
   
    def objective_fun(self, state):
        """The objective function (cost) at the current state"""
        self.update_rte_solver(state)
        self._loss = self.extinction_gradient_cost()[1]
        return self.loss
        
        
    def gradient(self, state):
        """The gradient at the current state"""
        self.update_rte_solver(state)
        parameter_gradient = [self.extinction_gradient_cost()[0]]
        state_gradient = np.concatenate(map(lambda grad: grad[self.gradient_mask], parameter_gradient))
        return state_gradient
    
    
    def save_ckpt(self):
        """Save a checkpoint of the Optimizer"""
        self._ckpt_time = time.time()
        timestr = time.strftime("%H%M%S")
        path = os.path.join(self.writer._tf_writer.log_dir,  timestr + '.ckpt')
        self.save(path)
            
            
    def callback(self, state):
        """
        The callback function invokes the callbacks defined by the writer (if any). 
        Additionally it keeps track of the iteration number 
        and if the checkpoint period has passed it saves a checkpoint of the Optimizer
        """
        self._iteration += 1
        
        # Writer callback functions
        if self.writer is not None:
            for fn in self.writer.callback_fns:
                fn()
        
        # Save checkpoint
        time_passed = time.time() - self._ckpt_time
        if (self.ckpt_period is not None) and (time_passed > self.ckpt_period):        
            self.save_ckpt()
        
        
    def init_minimization(self):
        """Initilize internal structures (masks) before starting the optimization process"""
        # Init masks for parameters
        self._num_parameters = []
        for param in self.parameters:
            if self.cloud_mask:
                param.set_mask(self.cloud_mask)
                self._num_parameters.append(param.num_parameters)
    
        # Init mask for gradient 
        if self.cloud_mask:
            if self.known_medium:
                grid = self.known_medium.grid + self.cloud_mask.grid
                self._gradient_mask = np.array(self.cloud_mask.resample(grid, method='nearest'), dtype=np.bool).ravel()
            else:
                self._gradient_mask = self.cloud_mask.data.ravel()
            

    def minimize(self, options, ckpt_period=None, method='L-BFGS-B'):
        """
        Minimize the cost function with respect to the parameters defined.
        
        Parameters
        ----------
        options: Dict
            The option dictionary
        ckpt_period: float
            Time in seconds between saving checkpoints. None means no checkpoints will be saved.
        method: str
            The optimization method.
        
        Notes
        -----
        Currently only L-BFGS-B optimization method is supported.
       
        For documentation: 
            https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html
        """
        
        self._ckpt_period = ckpt_period
        
        if method != 'L-BFGS-B':
            raise NotImplementedError('Optimization method not implemented')
        
        if self.iteration == 0:
            self.init_minimization() 

        initial_state = self.parameters_to_state()
        result = minimize(fun=self.objective_fun, 
                          x0=initial_state, 
                          method=method, 
                          jac=self.gradient,
                          bounds=self.get_bounds(),
                          options=options,
                          callback=self.callback)
        return result
    
    
    def parameters_to_state(self):
        """
        Transform parameters into a state by concatination
        
        Returns
        -------
        state: np.array()
            The state vector.
        """
        state = np.concatenate(map(lambda param: param.data[param.mask].ravel(), self.parameters))
        return state
    
    
    def get_bounds(self):
        """Transform shdom.Parameter bounds to a list of bounds that the scipy minimize expects."""
        bounds = []
        for param in self.parameters:
            bounds.extend(param.bounds)
        return bounds


    def save(self, path):
        """
        Save Optimizer to file.
        
        Parameters
        ----------
        path: str,
            Full path to file. 
        """
        params = self.__dict__.copy()
        params.pop('_writer')
        rte_solver = params.pop('_rte_solver')
        params['scene_params'] = rte_solver._scene_parameters
        params['numerical_params'] = rte_solver._numerical_parameters
        file = open(path, 'w')
        file.write(pickle.dumps(params, -1))
        file.close()


    def load(self, path):
        """
        Load Optimizer from file.
        
        Parameters
        ----------
        path: str,
            Full path to file. 
        """
        file = open(path, 'r')
        data = file.read()
        file.close()        
        params = pickle.loads(data)
        
        # Create an RTE solver object
        rte_solver = shdom.RteSolver(params.pop('scene_params'), 
                                     params.pop('numerical_params'))
        self.set_rte_solver(rte_solver)
        self.__dict__ = params
        
        
    @property
    def rte_solver(self):
        return self._rte_solver
    
    @property
    def radiances(self):
        return self._measurements.radiances
    
    @property
    def sensors(self):
        return self._measurements.sensors
    
    @property
    def parameters(self):
        return self._parameters
    
    @property
    def cloud_mask(self):
        return self._cloud_mask    
    
    @property
    def gradient_mask(self):
        return self._gradient_mask  
    
    @property
    def num_parameters(self):
        return self._num_parameters       
    
    @property
    def extinction_dependency(self):
        return self._extinction_dependency    
    
    @property
    def albedo_dependency(self):
        return self._albedo_dependency    

    @property
    def phase_dependency(self):
        return self._phase_dependency   
    
    @property
    def known_medium(self):
        return self._known_medium
    

    @property
    def writer(self):
        return self._writer      
    
    @property
    def iteration(self):
        return self._iteration
    
    @property
    def ckpt_period(self):
        return self._ckpt_period
    
    @property
    def loss(self):
        return self._loss        
    
    @property
    def images(self):
        return self._images      