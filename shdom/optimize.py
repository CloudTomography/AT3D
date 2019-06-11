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
from collections import OrderedDict


class GridPhaseEstimator(shdom.GridPhase):
    """TODO"""
    def __init__(self, legendre_table, index):
        super(GridPhaseEstimator, self).__init__(legendre_table, index)

   
class GridDataEstimator(shdom.GridData):
    """TODO"""
    def __init__(self, grid_data, min_bound=None, max_bound=None):
        super(GridDataEstimator, self).__init__(grid_data.grid, grid_data.data)
        self._min_bound = min_bound
        self._max_bound = max_bound
        self._mask = None
        self._num_parameters = self.init_num_parameters()
             
              
    def set_state(self, state):
        """TODO"""
        if self.mask is None:
            self._data = np.reshape(state, (self.shape))
        else:
            self._data[self.mask.data] = state
        
        
    def get_state(self):
        """TODO"""
        if self.mask is None:
            return self.data.ravel()
        else:
            return self.data[self.mask.data]
    
    
    def init_num_parameters(self):
        """TODO"""
        if self.mask is None:
            num_parameters = self.data.size
        else:
            num_parameters = np.count_nonzero(self.mask.data)
        return num_parameters


    def set_mask(self, mask):
        """
        TODO
        
        Parameters
        ----------
        mask: shdom.GridData object
            A boolean mask with True making cloudy voxels and False marking non-cloud region.     
        """
        self._mask = mask.resample(self.grid, method='nearest')
        self._data[np.bitwise_not(self._mask.data)] = 0.0
        self._num_parameters = self.init_num_parameters()
        
        
    def get_bounds(self):
        """TODO
        bounds to a list of bounds that the scipy minimize expects."""
        return [(self.min_bound, self.max_bound)] * self.num_parameters  
    
    
    def project_gradient(self, gradient):
        """TODO"""
        state_gradient = gradient.resample(self.grid)
        if self.mask is None:
            return state_gradient.data.ravel()
        else:
            return state_gradient.data[self.mask.data]
    
    
    @property
    def mask(self):
        return self._mask
    
    @property
    def num_parameters(self):
        return self._num_parameters
    
    @property
    def min_bound(self):
        return self._min_bound
    
    @property
    def max_bound(self):
        return self._max_bound
    
    
class ScattererEstimator(shdom.Scatterer):
    """TODO"""
    def __init__(self, extinction, albedo, phase):
        super(ScattererEstimator, self).__init__(extinction, albedo, phase)
        self._mask = None
        self._estimators = self.init_estimators()
        self._num_parameters = self.init_num_parameters()


    def init_estimators(self):
        """TODO"""
        return None
    
    
    def init_num_parameters(self):
        """TODO"""
        num_parameters = 0
        for estimator in self.estimators:
            num_parameters += estimator.num_parameters
        return num_parameters


    def set_mask(self, mask):
        """
        TODO
        
        Parameters
        ----------
        mask: shdom.GridData object
            A boolean mask with True making cloudy voxels and False marking non-cloud region.     
        """
        self._mask = mask.resample(self.grid, method='nearest')
        for estimator in self.estimators:
            estimator.set_mask(mask)
        self._num_parameters = self.init_num_parameters()
       
       
    def set_state(self, state):
        """TODO"""
        for estimator in self.estimators:
            estimator.set_state(state)        
        
        
    def get_state(self):
        """
        TODO
        """
        state = np.empty(shape=(0), dtype=np.float64)
        for estimator in self.estimators:
            state = np.concatenate((state, estimator.get_state()))
        return state


    def get_bounds(self):
        """TODO
        bounds to a list of bounds that the scipy minimize expects."""
        bounds = []
        for estimator in self.estimators:
            bounds.extend(estimator.get_bounds())
        return bounds
    
    
    def project_gradient(self, gradient):
        """TODO"""
        state_gradient = np.empty(shape=(0), dtype=np.float64)
        for estimator in self.estimators:
            state_gradient = np.concatenate((state_gradient, estimator.project_gradient(gradient)))
        return state_gradient
        
        
    @property
    def estimators(self):
        return self._estimators
    
    @property
    def num_parameters(self):
        return self._num_parameters

    @property
    def mask(self):
        return self._mask
    

class ScattererOpticalEstimator(shdom.ScattererEstimator):
    """TODO"""
    def __init__(self, extinction, albedo, phase):
        super(ScattererOpticalEstimator, self).__init__(extinction, albedo, phase)
        
    def init_estimators(self):
        """TODO"""
        estimators = []
        if self.extinction.__class__ == shdom.GridDataEstimator:
            estimators.append(self.extinction)
        if self.albedo.__class__ == shdom.GridDataEstimator:
            raise NotImplementedError("Albedo estimation not implemented")
        if self.phase.__class__ == shdom.GridPhaseEstimator:
            raise NotImplementedError("Phase estimation not implemented")           
        return estimators


class ScattererMicrophysicalEstimator(shdom.ScattererEstimator):
    """TODO"""
    def __init__(self, mie, lwc, reff, veff):
        self._mie = mie
        self._lwc = lwc
        self._reff = reff
        self._veff = veff
        super(ScattererMicrophysicalEstimator, self).__init__(
            extinction=mie.get_extinction(lwc, reff, veff), 
            albedo=mie.get_albedo(reff, veff), 
            phase=mie.get_phase(reff, veff, squeeze_table=True)
        )
    
    def init_estimators(self):
        """TODO"""
        estimators = []
        if self.lwc.__class__ == shdom.GridDataEstimator:
            estimators.append(self.lwc)
        if self.reff.__class__ == shdom.GridDataEstimator:
            estimators.append(self.reff)
        if self.veff.__class__ == shdom.GridDataEstimator:
            estimators.append(self.veff)   
        return estimators
    
    @property
    def mie(self):
        return self._mie
    
    @property
    def lwc(self):
        return self._lwc
    
    @property
    def reff(self):
        return self._reff   

    @property
    def veff(self):
        return self._veff  


class OpticalMediumEstimator(shdom.OpticalMedium):
    """TODO"""
    def __init__(self, grid=None):
        super(OpticalMediumEstimator, self).__init__(grid)
        self._estimators = OrderedDict()
        self._num_parameters = []
        self._num_unknown_scatterers = 0
        self._unknown_scatterers_indices = []
        
    def add_scatterer(self, scatterer, name=None):
        """TODO"""
        super(OpticalMediumEstimator, self).add_scatterer(scatterer, name)
        if scatterer.__class__ == shdom.ScattererEstimator:
            name = 'scatterer{:d}'.format(self._num_scatterers) if name is None else name 
            self._estimators[name] = scatterer
            self._num_parameters.append(scatterer.num_parameters)
            self._unknown_scatterers_indices.append(self.num_scatterers)
            self._num_unknown_scatterers += 1
            
            
    def set_state(self, state):
        """TODO"""
        states = np.split(state, np.cumsum(self.num_parameters[:-1]))
        for (name, estimator), state in zip(self.estimators.iteritems(), states):
            estimator.set_state(state)
            self.update_scatterer(name, estimator)
    
    
    def get_state(self):
        """TODO"""
        state = np.empty(shape=(0),dtype=np.float64)
        for estimator in self.estimators.itervalues():
            state = np.concatenate((state, estimator.get_state()))
        return state


    def get_bounds(self):
        """TODO
        bounds to a list of bounds that the scipy minimize expects."""
        bounds = []
        for estimator in self.estimators.itervalues():
            bounds.extend(estimator.get_bounds())
        return bounds
    
    
    def compute_gradient(self, rte_solver, projection, measurements):
        """The objective function (cost) at the current state"""
        
        multiview = projection.__class__ is shdom.sensor.MultiViewProjection
        
        if isinstance(projection.npix, list):
            total_pix = np.sum(projection.npix)
        else:
            total_pix = projection.npix
    
        gradient, loss, radiance = core.gradient(
            partder=self.unknown_scatterers_indices,
            numder=self.num_unknown_scatterers,
            extflag=self.extflag,
            phaseflag=self.phaseflag,
            phaseder=np.zeros_like(rte_solver._legen),             
            nstokes=rte_solver._nstokes,
            nstleg=rte_solver._nstleg,
            nx=rte_solver._nx,
            ny=rte_solver._ny,
            nz=rte_solver._nz,
            bcflag=rte_solver._bcflag,
            ipflag=rte_solver._ipflag,   
            npts=rte_solver._npts,
            nbpts=rte_solver._nbpts,
            ncells=rte_solver._ncells,
            nbcells=rte_solver._nbcells,
            ml=rte_solver._ml,
            mm=rte_solver._mm,
            ncs=rte_solver._ncs,
            nlm=rte_solver._nlm,
            numphase=rte_solver._pa.numphase,
            nmu=rte_solver._nmu,
            nphi0max=rte_solver._nphi0max,
            nphi0=rte_solver._nphi0,
            maxnbc=rte_solver._maxnbc,
            ntoppts=rte_solver._ntoppts,
            nbotpts=rte_solver._nbotpts,
            nsfcpar=rte_solver._nsfcpar,
            gridptr=rte_solver._gridptr,
            neighptr=rte_solver._neighptr,
            treeptr=rte_solver._treeptr,             
            shptr=rte_solver._shptr,
            bcptr=rte_solver._bcptr,
            cellflags=rte_solver._cellflags,
            iphase=rte_solver._iphase[:rte_solver._npts],
            deltam=rte_solver._deltam,
            solarmu=rte_solver._solarmu,
            solaraz=rte_solver._solaraz,
            gndtemp=rte_solver._gndtemp,
            gndalbedo=rte_solver._gndalbedo,
            skyrad=rte_solver._skyrad,
            waveno=rte_solver._waveno,
            wavelen=rte_solver._wavelen,
            mu=rte_solver._mu,
            phi=rte_solver._phi.reshape(rte_solver._nmu, -1),
            wtdo=rte_solver._wtdo.reshape(rte_solver._nmu, -1),
            xgrid=rte_solver._xgrid,
            ygrid=rte_solver._ygrid,
            zgrid=rte_solver._zgrid,
            gridpos=rte_solver._gridpos,
            sfcgridparms=rte_solver._sfcgridparms,
            bcrad=rte_solver._bcrad,
            extinct=rte_solver._extinct[:rte_solver._npts],
            albedo=rte_solver._albedo[:rte_solver._npts],
            legen=rte_solver._legen.reshape(rte_solver._nleg+1, -1),            
            dirflux=rte_solver._dirflux,
            fluxes=rte_solver._fluxes,
            source=rte_solver._source,
            camx=projection.x,
            camy=projection.y,
            camz=projection.z,
            cammu=projection.mu,
            camphi=projection.phi,
            npix=total_pix,       
            srctype=rte_solver._srctype,
            sfctype=rte_solver._sfctype,
            units=rte_solver._units,
            measurements=measurements,
            rshptr=rte_solver._rshptr,
            radiance=rte_solver._radiance,
            total_ext=rte_solver._total_ext[:rte_solver._npts])       
        
        # Project rte_solver base grid gradient to the state space
        gradient = gradient.reshape(self.grid.shape + tuple([self.num_unknown_scatterers]))
        state_gradient = np.empty(shape=(0), dtype=np.float64)
        for i, estimator in enumerate(self.estimators.values()):
            state_gradient = np.concatenate(
                (state_gradient, estimator.project_gradient(GridData(self.grid, gradient[...,i])))
            )

        # Split images into different sensor images
        if multiview:
            split_indices = np.cumsum(projection.npix[:-1])        
            images = np.split(radiance, split_indices)
            images = [
                image.reshape(resolution, order='F') 
                for image, resolution in zip(images, projection.resolution) 
            ]  
            
        return state_gradient, loss, images

    @property
    def estimators(self):
        return self._estimators
    
    @property
    def num_parameters(self):
        return self._num_parameters
        
    @property
    def num_unknown_scatterers(self):
        return self._num_unknown_scatterers
    
    @property
    def unknown_scatterers_indices(self):
        return np.array(self._unknown_scatterers_indices, dtype=np.int32)
    
    
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
        
        
    def monitor_parameter_error(self, ground_truth, ckpt_period=-1):
        """
        Monitor relative and overall mass error (epsilon, delta) as defined by:
          Amit Aides et al, "Multi sky-view 3D aerosol distribution recovery".
        
        Parameters
        ----------
        ground_truth: GridData
            A GridData of the ground truth parameters.
        ckpt_period: float
           time [seconds] between updates. setting ckpt_period=-1 will log at every iteration.
        """
        self._param_ckpt_period = ckpt_period
        self._param_ckpt_time = time.time()        
        self._ground_truth_params = ground_truth
        self._callback_fns.append(self.param_error)
        
        
    def param_error(self):
        """Callback function the is called every optimizer iteration parameter error monitoring is set."""
        time_passed = time.time() - self._image_ckpt_time 
        if time_passed > self._param_ckpt_period:
            est_param = self.optimizer.medium.get_scatterer('cloud estimator').extinction
            gt_param = self._ground_truth_params
            delta = (np.linalg.norm(est_param.data.ravel(),1) - np.linalg.norm(gt_param.data.ravel(),1)) / np.linalg.norm(gt_param.data.ravel(),1)
            epsilon = np.linalg.norm((est_param - gt_param).data.ravel(), 1) / np.linalg.norm(gt_param.data.ravel(),1)
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
        
        if measurements.camera.projection.__class__ == shdom.MultiViewProjection:
            self._projections = measurements.camera.projection.projection_list
        else:
            self._projections = [measurements.camera.projection]
        self._images = measurements.images
        
        
    def carve(self, grid, thresholds, agreement=0.75):
        """
        Carves out the cloud geometry on the grid. 
        A threshold on the radiances is found by a cv2 adaptive image threshold.
        
        Parameters
        ----------
        grid: shdom.Grid
            A grid object.
        thresholds: list or float
            Either a constant threshold or a list of len(thresholds)=num_projections is used as for masking.
        agreement: float
            the precentage of pixels that should agree on a cloudy voxels to set it to True in the mask
        
        Returns
        -------
        mask: shdom.GridData object
            A boolean mask with True marking cloudy voxels and False marking non-cloud region.
        """
        
        self._rte_solver.set_grid(grid)
        volume = np.zeros((grid.nx, grid.ny, grid.nz))
        
        thresholds = np.array(thresholds)
        if thresholds.size == 1:
            thresholds = np.repeat(thresholds, len(self._images))
        else:
            assert thresholds.size == len(self._images), 'thresholds (len={}) should be of the same' \
                   'length as the number of images (len={})'.format(thresholds.size,  len(self._images))
            
        for projection, image, threshold in zip(self._projections, self._images, thresholds):
            
            image_mask = image > threshold
                
            projection = projection[image_mask.ravel(order='F') == 1]
            
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
                camx=projection.x,
                camy=projection.y,
                camz=projection.z,
                cammu=projection.mu,
                camphi=projection.phi,
                npix=projection.npix,
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
       [required] optimizer.set_medium() 
       [optional] optimizer.set_writer() 
    
    Notes
    -----
    Currently only Extinction optimization is supported.
    """
    def __init__(self):
        self._medium = None
        self._rte_solver = None
        self._measurements = None
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
    
    
    def set_medium_estimator(self, medium):
        """
        Set the MediumEstimator for the optimizer.
        
        Parameters
        ----------
        medium: shdom.MediumEstimator
            The MediumEstimator
        """
        self._medium = medium       


    def set_rte_solver(self, rte_solver):
        """
        Set the RteSolver for the SHDOM iterations.
        
        Parameters
        ----------
        rte_solver: shdom.RteSolver
            The RteSolver
        """
        self._rte_solver = rte_solver


    def set_writer(self, writer):
        """
        Set a log writer to upload summaries into tensorboard.
        
        Parameters
        ----------
        writer: shdom.SummaryWriter
            Wrapper for the tensorboardX summary writer.
        """
        self._writer = writer
        if writer is not None:
            self._writer.attach_optimizer(self)
        
   
    def objective_fun(self, state):
        """The objective function (cost) at the current state"""
        self.set_state(state)
        self.rte_solver.solve(maxiter=100, verbose=False)
        gradient, loss, images = self.medium.compute_gradient(rte_solver=self.rte_solver, 
                                                              projection=self.projection, 
                                                              measurements=self.radiances)
        self._loss = loss
        self._images = images
        return loss, gradient
    
    
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
            self.rte_solver.init_medium(self.medium)
            self._num_parameters = self.medium.num_parameters
            
        result = minimize(fun=self.objective_fun, 
                          x0=self.get_state(), 
                          method=method, 
                          jac=True,
                          bounds=self.get_bounds(),
                          options=options,
                          callback=self.callback)
        return result
    
    
    def get_bounds(self):
        """TODO"""
        return self.medium.get_bounds()
    
    
    def get_state(self):
        """TODO"""
        return self.medium.get_state()
    
    
    def set_state(self, state):
        """TODO"""
        self.medium.set_state(state)
        self.rte_solver.init_medium(self.medium)
        
        
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
    def medium(self):
        return self._medium    
    
    @property
    def radiances(self):
        return self._measurements.radiances
    
    @property
    def projection(self):
        return self._measurements.camera.projection
    
    
    @property
    def num_parameters(self):
        return self._num_parameters

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