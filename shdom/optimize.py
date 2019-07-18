"""
Optimization and related objects to monitor and log the optimization proccess.
"""

import shdom
import numpy as np
import core, time, os, copy
from scipy.optimize import minimize
from shdom import GridData
import dill as pickle
import cv2
import itertools
from joblib import Parallel, delayed
from collections import OrderedDict
from scipy.interpolate import RegularGridInterpolator


class SizeDistributionDerivative(shdom.SizeDistribution):
    """TODO"""
    def __init__(self, type='gamma'):
        super(SizeDistributionDerivative, self).__init__(type)
        
    def compute_nd(self, type, radii, particle_density=1.0):
        """TODO"""
        self._radii = radii
        self._nsize = radii.shape[0]           
        self._pardens = particle_density        
          
        if type=='reff':
            self._nd = self._compute_re_derivative()
        elif type=='veff':
            self._nd = self._compute_ve_derivative()
        else:
            raise AttributeError('Derivative type is either reff or veff')
        
        nd = self.nd.T.reshape((self.nretab, self.nvetab, self.nsize), order='F')
        self._nd_interpolator = RegularGridInterpolator(
            (self.reff, self.veff), nd, bounds_error=False, fill_value=0.0)
        
    def _compute_re_derivative(self):
        """TODO"""
        reff, alpha = np.meshgrid(self.reff, self.alpha)  
        eps = 1e-5
        nd1 = core.make_multi_size_dist(
            distflag=self.distflag,
            pardens=self.pardens,
            nsize=self.nsize,
            radii=self.radii,
            reff=reff.ravel(),
            alpha=alpha.ravel(),
            gamma=self.gamma,
            ndist=reff.size)
        nd2 = core.make_multi_size_dist(
            distflag=self.distflag,
            pardens=self.pardens,
            nsize=self.nsize,
            radii=self.radii,
            reff=reff.ravel() + eps,
            alpha=alpha.ravel(),
            gamma=self.gamma,
            ndist=reff.size)
        return (nd2 - nd1) / eps
        
    def _compute_ve_derivative(self):
        """TODO"""
        reff, veff = np.meshgrid(self.reff, self.veff)  
        eps = 1e-5
        
        if self.distflag == 'G':
            alpha1 = 1.0 / veff.ravel() - 3.0
            alpha2 = 1.0 / (veff.ravel() + eps) - 3.0
        if self.distflag == 'L':
            alpha1 = np.sqrt(np.log(veff.ravel() + 1.0))
            alpha2 = np.sqrt(np.log(veff.ravel() + eps + 1.0))
            
        nd1 = core.make_multi_size_dist(
            distflag=self.distflag,
            pardens=self.pardens,
            nsize=self.nsize,
            radii=self.radii,
            reff=reff.ravel(),
            alpha=alpha1,
            gamma=self.gamma,
            ndist=reff.size)
        nd2 = core.make_multi_size_dist(
            distflag=self.distflag,
            pardens=self.pardens,
            nsize=self.nsize,
            radii=self.radii,
            reff=reff.ravel(),
            alpha=alpha2,
            gamma=self.gamma,
            ndist=reff.size)
        return (nd2 - nd1) / eps
    
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
    
    
class ScattererEstimator(object):
    """TODO"""
    def __init__(self):      
        self._mask = None         
        self._estimators = self.init_estimators()
        self._derivatives = self.init_derivatives()
        self._num_parameters = self.init_num_parameters()
        
    def init_estimators(self):
        """TODO"""
        return OrderedDict()
    
    def init_derivatives(self):
        """TODO"""
        return OrderedDict()    

    def init_num_parameters(self):
        """TODO"""
        num_parameters = []
        for estimator in self.estimators.itervalues():
            num_parameters.append(estimator.num_parameters)
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
        for estimator in self.estimators.itervalues():
            estimator.set_mask(mask)
        self._num_parameters = self.init_num_parameters()
       
    def set_state(self, state):
        """TODO"""
        states = np.split(state, np.cumsum(self.num_parameters[:-1]))
        for estimator, state in zip(self.estimators.itervalues(), states):
            estimator.set_state(state)        
    
    def get_state(self):
        """
        TODO
        """
        state = np.empty(shape=(0), dtype=np.float64)
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
        
        
    def project_gradient(self, gradient):
        """TODO"""
        state_gradient = np.empty(shape=(0), dtype=np.float64)
        for estimator in self.estimators.itervalues():
            state_gradient = np.concatenate((state_gradient, estimator.project_gradient(gradient)))
        return state_gradient
        
        
    @property
    def estimators(self):
        return self._estimators
    
    @property
    def derivatives(self):
        return self._derivatives    
    
    @property
    def num_parameters(self):
        return self._num_parameters

    @property
    def mask(self):
        return self._mask

    
    
class OpticalScattererEstimator(shdom.OpticalScatterer, ScattererEstimator):
    """TODO"""
    def __init__(self, wavelength, extinction, albedo, phase):
        shdom.OpticalScatterer.__init__(self, wavelength, extinction, albedo, phase)
        ScattererEstimator.__init__(self)
        
    def init_estimators(self):
        """TODO"""
        estimators = OrderedDict()
        if isinstance(self.extinction, shdom.GridDataEstimator):
            estimators['extinction'] = self.extinction
        if isinstance(self.albedo, shdom.GridDataEstimator):
            raise NotImplementedError("Albedo estimation not implemented")
        if isinstance(self.phase, shdom.GridPhaseEstimator):
            raise NotImplementedError("Phase estimation not implemented")           
        return estimators

    def init_derivatives(self):
        """TODO"""
        derivatives = OrderedDict()
        if isinstance(self.extinction, shdom.GridDataEstimator):
            extinction = shdom.GridData(self.extinction.grid, np.ones_like(self.extinction.data))
            albedo = shdom.GridData(self.albedo.grid, np.zeros_like(self.albedo.data))
            legen_table = shdom.LegendreTable(np.zeros((self.phase.legendre_table.maxleg+1), dtype=np.float32), table_type=self.phase.legendre_table.table_type)
            phase = shdom.GridPhase(legen_table, shdom.GridData(self.phase.index.grid, np.ones_like(self.phase.index.data)))
            derivatives['extinction'] = shdom.OpticalScatterer(self.wavelength, extinction, albedo, phase)
        if isinstance(self.albedo, shdom.GridDataEstimator):
            raise NotImplementedError("Albedo estimation not implemented")
        if isinstance(self.phase, shdom.GridPhaseEstimator):
            raise NotImplementedError("Phase estimation not implemented")           
        return derivatives


class MicrophysicalScattererEstimator(shdom.MicrophysicalScatterer, ScattererEstimator):
    """TODO"""
    def __init__(self, lwc, reff, veff):
        shdom.MicrophysicalScatterer.__init__(self, lwc, reff, veff)
        ScattererEstimator.__init__(self)
        self._re_derivative = OrderedDict()
        self._ve_derivative = OrderedDict()
        
    def init_estimators(self):
        """TODO"""
        estimators = OrderedDict()
        if isinstance(self.lwc, shdom.GridDataEstimator):
            estimators['lwc'] = self.lwc
        if isinstance(self.reff, shdom.GridDataEstimator):
            estimators['reff'] = self.reff
        if isinstance(self.veff, shdom.GridDataEstimator):
            estimators['veff'] = self.veff
        return estimators
    
    
    def add_mie(self, mie):
        """TODO"""
        super(MicrophysicalScattererEstimator, self).add_mie(mie)
        if self.estimators.has_key('reff'):
            self.add_re_derivative(mie)
        if self.estimators.has_key('veff'):
            self.add_ve_derivative(mie)
        
        
    def add_re_derivative(self, mie):
        """TODO"""
        if isinstance(mie, shdom.MiePolydisperse):
            mie_list = [mie]
        elif isinstance(mie, dict):
            mie_list = mie.values()
        
        for mie in mie_list:
            self._re_derivative[mie.wavelength] = mie.get_re_derivative()


    def add_ve_derivative(self, mie):
        """TODO"""
        raise NotImplementedError
    
        
class MediumEstimator(shdom.Medium):
    """TODO"""
    def __init__(self, grid=None):
        super(MediumEstimator, self).__init__(grid)
        self._estimators = OrderedDict()
        self._num_parameters = []
        self._unknown_scatterers_indices =  np.empty(shape=(0),dtype=np.int32)
        self._num_derivatives = 0
        
    def add_scatterer(self, scatterer, name=None):
        """TODO"""
        super(MediumEstimator, self).add_scatterer(scatterer, name)
        if issubclass(type(scatterer), shdom.ScattererEstimator):
            name = 'scatterer{:d}'.format(self._num_scatterers) if name is None else name 
            num_estimators = len(scatterer.estimators)
            self._estimators[name] = scatterer
            self._num_parameters.append(np.sum(scatterer.num_parameters))
            self._unknown_scatterers_indices = np.concatenate((
                self.unknown_scatterers_indices, 
                np.full(num_estimators, self.num_scatterers, dtype=np.int32)))
            self._num_derivatives += num_estimators
                 
    def set_state(self, state):
        """TODO"""
        states = np.split(state, np.cumsum(self.num_parameters[:-1]))
        for (name, estimator), state in zip(self.estimators.iteritems(), states):
            estimator.set_state(state)
            self.scatterers[name] = estimator
    
    
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
    
    
    def get_derivatives(self, rte_solver):
        """TODO"""
        dext = np.zeros(shape=[rte_solver._nbpts, self.num_derivatives], dtype=np.float32)
        dalb = np.zeros(shape=[rte_solver._nbpts, self.num_derivatives], dtype=np.float32)
        diphase = np.zeros(shape=[rte_solver._nbpts, self.num_derivatives], dtype=np.float32)
    
        i=0
        for estimator in self.estimators.itervalues():
            for derivative in estimator.derivatives.itervalues():
                if isinstance(derivative, shdom.MicrophysicalScatterer) or isinstance(derivative, shdom.MultispectralScatterer):
                    derivative = derivative.get_optical_scatterer(rte_solver._wavelen)            
                resampled_derivative = derivative.resample(self.grid)
                extinction = resampled_derivative.extinction.data
                albedo = resampled_derivative.albedo.data
                iphase = resampled_derivative.phase.iphasep                
                if rte_solver._bcflag == 0:
                    extinction = np.pad(extinction, ((0,1),(0,1),(0,0)), 'wrap')
                    albedo = np.pad(albedo, ((0,1),(0,1),(0,0)), 'wrap')
                    iphase = np.pad(iphase, ((0,1),(0,1),(0,0)), 'wrap')
                elif rte_solver._bcflag == 1:
                    extinction = np.pad(extinction, ((0,0),(0,1),(0,0)), 'wrap')
                    albedo = np.pad(albedo, ((0,0),(0,1),(0,0)), 'wrap')
                    iphase = np.pad(iphase, ((0,0),(0,1),(0,0)), 'wrap')
                elif rte_solver._bcflag == 2:
                    extinction = np.pad(extinction, ((0,1),(0,0),(0,0)), 'wrap')
                    albedo = np.pad(albedo, ((0,1),(0,0),(0,0)), 'wrap')
                    iphase = np.pad(iphase, ((0,1),(0,0),(0,0)), 'wrap')            
    
                dext[:,i] = extinction.ravel()
                dalb[:,i] = albedo.ravel()
                diphase[:,i] = iphase.ravel()      
    
                if i == 0:
                    leg_table = copy.deepcopy(resampled_derivative.phase.legendre_table)
                else:
                    leg_table.append(copy.deepcopy(resampled_derivative.phase.legendre_table))                
                i += 1
                
        dleg = leg_table.data  
        dnumphase= leg_table.numphase
        dphasetab = core.precompute_phase_check(
            nscatangle=rte_solver._nscatangle,
            numphase=dnumphase,
            ml=rte_solver._ml,
            nlm=rte_solver._nlm,
            nleg=rte_solver._nleg,
            legen=dleg,
            deltam=False
        )                
        return dext, dalb, diphase, dleg, dphasetab, dnumphase
                 
                    
    def core_grad(self, rte_solver, projection, radiance):
        """
        TODO
        """
        
        if isinstance(projection.npix, list):
            total_pix = np.sum(projection.npix)
        else:
            total_pix = projection.npix

        gradient, loss, radiance = core.gradient(
            partder=self.unknown_scatterers_indices,
            numder=self.num_derivatives,
            dext=rte_solver._dext,
            dalb=rte_solver._dalb,
            diphase=rte_solver._diphase,
            dleg=rte_solver._dleg,
            dphasetab=rte_solver._dphasetab,
            dnumphase=rte_solver._dnumphase,
            nscatangle=rte_solver._nscatangle,
            phasetab=rte_solver._phasetab,
            ylmsun=rte_solver._ylmsun,
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
            measurements=radiance,
            rshptr=rte_solver._rshptr,
            radiance=rte_solver._radiance,
            total_ext=rte_solver._total_ext[:rte_solver._npts]
        )
        
        return gradient, loss, radiance
        
    def compute_gradient(self, rte_solver, measurements, n_jobs):
        """
        The objective function (cost) and gradient at the current state.
        
        TODO 
        
        Parameters
        ----------
        n_jobs: int,
            The number of jobs to divide the gradient computation into.     
        """
                
        # If rendering several atmospheres (e.g. multi-spectral rendering)
        if isinstance(rte_solver, shdom.RteSolverArray):
            num_channels = rte_solver.num_solvers
            rte_solvers = rte_solver
        else:
            num_channels = 1
            rte_solvers = [rte_solver]
        
        # Pre-computation of phase-function and derivatives for all solvers.
        for rte_solver in rte_solvers:
            rte_solver._phasetab = core.precompute_phase_check(
                nscatangle=rte_solver._nscatangle,
                numphase=rte_solver._pa.numphase,
                ml=rte_solver._ml,
                nlm=rte_solver._nlm,
                nleg=rte_solver._nleg,
                legen=rte_solver._legen.reshape(rte_solver._nleg+1, -1),
                deltam=rte_solver._deltam
            )
            rte_solver._dext, rte_solver._dalb, rte_solver._diphase, \
                rte_solver._dleg, rte_solver._dphasetab, rte_solver._dnumphase = self.get_derivatives(rte_solver)
            
            for estimator in self.estimators.itervalues():
                for derivative in estimator.derivatives.itervalues():
                    if isinstance(derivative, shdom.MicrophysicalScatterer) or isinstance(derivative, shdom.MultispectralScatterer):
                        derivative = derivative.get_optical_scatterer(wavelength)                    
                        
                        
                
        projection = measurements.camera.projection
        sensor = measurements.camera.sensor
        radiances = measurements.radiances
        
        # Sequential or parallel processing using multithreading (threadsafe Fortran)
        if n_jobs > 1:           
            output = Parallel(n_jobs=n_jobs, backend="threading", verbose=0)(
                delayed(self.core_grad, check_pickle=False)(
                    rte_solver=rte_solvers[channel],
                    projection=projection,
                    radiance=spectral_radiance[..., channel]) for channel, (projection, spectral_radiance) in 
                itertools.product(range(num_channels), zip(projection.split(n_jobs), np.array_split(radiances, n_jobs)))
            )
        else:
            output = [self.core_grad(rte_solvers[channel], projection, radiances[...,channel]) for channel in range(num_channels)]
        
        # Sum over all the losses of the different channels
        loss = np.sum(map(lambda x: x[1], output))
        
        # Project rte_solver base grid gradient to the state space
        gradient = sum(map(lambda x: x[0], output))
        gradient = gradient.reshape(self.grid.shape + tuple([self.num_derivatives]))
        images = sensor.make_images(np.concatenate(map(lambda x: x[2], output)), projection, num_channels)
        
        state_gradient = np.empty(shape=(0), dtype=np.float64)
        for i, estimator in enumerate(self.estimators.values()):
            state_gradient = np.concatenate(
                (state_gradient, estimator.project_gradient(GridData(self.grid, gradient[...,i])))
            )
            
        images = sensor.make_images(np.concatenate(map(lambda x: x[2], output)), projection, num_channels)
        
        return state_gradient, loss, images


    @property
    def estimators(self):
        return self._estimators
    
    @property
    def num_parameters(self):
        return self._num_parameters
    
    @property
    def num_parameters(self):
        return self._num_parameters    

    @property
    def unknown_scatterers_indices(self):
        return self._unknown_scatterers_indices

    @property
    def num_derivatives(self):
        return self._num_derivatives 
    
    
    
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
        i = 0
        for view, image in enumerate(acquired_images):
            
            # for multispectral images
            if image.ndim == 3:
                for channel in range(image.shape[2]):
                    self._image_vmax.append(image[...,channel].max() * 1.25)
                    self.tf_writer.add_image(
                        tag='Acquiered Image view {} channel {}'.format(view, channel), 
                        img_tensor=(image[...,channel] / self._image_vmax[i]),
                        dataformats='HW'
                    )
                    i += 1
            # for single wavelength       
            else:
                self._image_vmax.append(image.max() * 1.25)
                self.tf_writer.add_image(
                    tag='Acquiered Image view {}'.format(view), 
                    img_tensor=(image / self._image_vmax[i]),
                    dataformats='HW'
                )
                i += 1
        self._callback_fns.append(self.images)
    
    
    def images(self):
        """Callback function the is called every optimizer iteration image monitoring is set."""
        time_passed = time.time() - self._image_ckpt_time 
        if time_passed > self._image_ckpt_period:
            self._image_ckpt_time  = time.time()
            
            i = 0
            for view, image in enumerate(self.optimizer.images):
                
                # for multispectral images
                if image.ndim == 3:
                    for channel in range(image.shape[2]):
                        self.tf_writer.add_image(
                            tag='Retrieval Image view {} channel {}'.format(view, channel), 
                            img_tensor=(image[...,channel] / self._image_vmax[i]),
                            dataformats='HW'
                        )
                        i +=1
                # for single wavelength    
                else:
                    self.tf_writer.add_image(
                        tag='Retrieval Image view {}'.format(view), 
                        img_tensor=(image / self._image_vmax[i]),
                        dataformats='HW'
                    ) 
                    i += 1


    def monitor_loss(self):
        """Monitor the loss at every iteration."""
        self._callback_fns.append(self.loss)
    
    
    def loss(self):
        """Callback function the is called every optimizer iteration loss monitoring is set."""
        self.tf_writer.add_scalar('loss', self.optimizer.loss, self.optimizer.iteration)
        
        
    def monitor_scatterer_error(self, estimated_scatterer_name, ground_truth_scatterer, ckpt_period=-1):
        """
        Monitor relative and overall mass error (epsilon, delta) as defined at:
          Amit Aides et al, "Multi sky-view 3D aerosol distribution recovery".
        
        Parameters
        ----------
        estimated_scatterer_name: str
            The name of the scatterer to monitor
        ground_truth_scatterer: shdom.Scatterer
            The ground truth medium.
        ckpt_period: float
           time [seconds] between updates. setting ckpt_period=-1 will log at every iteration.
        """
        self._param_ckpt_period = ckpt_period
        self._param_ckpt_time = time.time()
        if hasattr(self, '_ground_truth_scatterer'):
            self._ground_truth_scatterer[estimated_scatterer_name] = ground_truth_scatterer
        else:
            self._ground_truth_scatterer = OrderedDict({estimated_scatterer_name: ground_truth_scatterer})
        self._callback_fns.append(self.param_error)
        
        
    def param_error(self):
        """Callback function the is called every optimizer iteration parameter error monitoring is set."""
        time_passed = time.time() - self._image_ckpt_time 
        if time_passed > self._param_ckpt_period:
            for scatterer_name, gt_scatterer in self._ground_truth_scatterer.iteritems():
                est_scatterer = self.optimizer.medium.get_scatterer(scatterer_name)
                for parameter_name, parameter in est_scatterer.estimators.iteritems():
                    est_param = parameter.data[parameter.mask.data].ravel()
                    gt_param = getattr(gt_scatterer, parameter_name).data[parameter.mask.data].ravel()
                    delta = (np.linalg.norm(est_param, 1) - np.linalg.norm(gt_param, 1)) / np.linalg.norm(gt_param, 1)
                    epsilon = np.linalg.norm((est_param - gt_param), 1) / np.linalg.norm(gt_param,1)
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
        
        if isinstance(measurements.camera.projection, shdom.MultiViewProjection):
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
       [required] optimizer.set_medium_estimator() 
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
    
    
    def set_medium_estimator(self, medium_estimator):
        """
        Set the MediumEstimator for the optimizer.
        
        Parameters
        ----------
        medium_estimator: shdom.MediumEstimator
            The MediumEstimator
        """
        self._medium = medium_estimator       


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
        gradient, loss, images = self.medium.compute_gradient(
            rte_solver=self.rte_solver,
            measurements=self.measurements,
            n_jobs=self._n_jobs
        )
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
        

    def minimize(self, options, ckpt_period=None, method='L-BFGS-B', n_jobs=1):
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
        n_jobs: int, default=1
            The number of jobs to divide the gradient computation into.
        
        Notes
        -----
        Currently only L-BFGS-B optimization method is supported.
       
        For documentation: 
            https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html
        """
        
        self._ckpt_period = ckpt_period
        self._n_jobs = n_jobs
        if method != 'L-BFGS-B':
            raise NotImplementedError('Optimization method not implemented')
        
        if self.iteration == 0:
            self.init_optimizer()
            
        result = minimize(fun=self.objective_fun, 
                          x0=self.get_state(), 
                          method=method, 
                          jac=True,
                          bounds=self.get_bounds(),
                          options=options,
                          callback=self.callback)
        return result
    
    
    def init_optimizer(self):
        """TODO"""
        self.rte_solver.init_medium(self.medium)
        self._num_parameters = self.medium.num_parameters   
        
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
    def measurements(self):
        return self._measurements
    
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