"""
Optimization and related objects to monitor and log the optimization proccess.
"""
import numpy as np
import time, os, copy
from scipy.optimize import minimize
import shdom
from shdom import GridData, core, float_round
import dill as pickle
import itertools
from joblib import Parallel, delayed
from collections import OrderedDict


class OpticalScattererDerivative(shdom.OpticalScatterer):
    """TODO"""
    def __init__(self, wavelength, extinction=None, albedo=None, phase=None):
        super(OpticalScattererDerivative, self).__init__(wavelength, extinction, albedo, phase)
        
    def resample(self, grid):
        """TODO"""
        extinction = self.extinction.resample(grid)
        albedo = self.albedo.resample(grid)
        phase = self.phase.resample(grid)            
        return shdom.OpticalScattererDerivative(self.wavelength, extinction, albedo, phase)
    
    @property
    def extinction(self):
        return self._extinction
    
    @extinction.setter
    def extinction(self, val):
        self._extinction = val
    
    @property
    def albedo(self):
        return self._albedo    
    
    @albedo.setter
    def albedo(self, val):
        self._albedo = val
        
    
        
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
        for estimator in self.estimators.values():
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
        for estimator in self.estimators.values():
            estimator.set_mask(mask)
        self._num_parameters = self.init_num_parameters()
       
    def set_state(self, state):
        """TODO"""
        states = np.split(state, np.cumsum(self.num_parameters[:-1]))
        for estimator, state in zip(self.estimators.values(), states):
            estimator.set_state(state)        
    
    def get_state(self):
        """
        TODO
        """
        state = np.empty(shape=(0), dtype=np.float64)
        for estimator in self.estimators.values():
            state = np.concatenate((state, estimator.get_state()))
        return state


    def get_bounds(self):
        """TODO
        bounds to a list of bounds that the scipy minimize expects."""
        bounds = []
        for estimator in self.estimators.values():
            bounds.extend(estimator.get_bounds())
        return bounds
        
        
    def project_gradient(self, gradient):
        """TODO"""
        state_gradient = np.empty(shape=(0), dtype=np.float64)
        for estimator in self.estimators.values():
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
            derivatives['extinction'] = self.init_extinction_derivative()
        if isinstance(self.albedo, shdom.GridDataEstimator):
            derivatives['albedo'] = self.init_albedo_derivative()
        if isinstance(self.phase, shdom.GridPhaseEstimator):
            derivatives['phase'] = self.init_phase_derivative()   
        return derivatives

    def init_extinction_derivative(self):
        """TODO"""
        extinction = shdom.GridData(self.extinction.grid, np.ones_like(self.extinction.data))
        albedo = shdom.GridData(self.albedo.grid, np.zeros_like(self.albedo.data))
        legen_table = shdom.LegendreTable(np.zeros((self.phase.legendre_table.maxleg+1), dtype=np.float32), table_type=self.phase.legendre_table.table_type)
        phase = shdom.GridPhase(legen_table, shdom.GridData(self.phase.index.grid, np.ones_like(self.phase.index.data)))
        derivative = shdom.OpticalScattererDerivative(self.wavelength, extinction, albedo, phase)        
        return derivative
        
    def init_albedo_derivative(self):
        """TODO"""
        raise NotImplementedError("Albedo estimation not implemented")
    
    
    def init_phase_derivative(self):
        """TODO"""
        raise NotImplementedError("Phase estimation not implemented")  

    def get_derivative(self, derivative_type, wavelength):
        """TODO"""
        if derivative_type == 'extinction':
            derivative = self.derivatives['extinction']
        elif derivative_type == 'albedo':
            derivative = self.derivatives['albedo']
        elif derivative_type == 'phase':
            derivative = self.derivatives['phase']   
        else:
            raise AttributeError('derivative type {} not supported'.format(derivative_type))
        return derivative
        
        
class MicrophysicalScattererEstimator(shdom.MicrophysicalScatterer, ScattererEstimator):
    """TODO"""
    def __init__(self, mie, lwc, reff, veff):
        shdom.MicrophysicalScatterer.__init__(self, lwc, reff, veff)
        self.add_mie(mie)
        ScattererEstimator.__init__(self)
        
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
    
    
    def init_derivatives(self):
        """TODO"""
        derivatives = OrderedDict()
        if isinstance(self.lwc, shdom.GridDataEstimator):
            derivatives['lwc'] = self.init_lwc_derivative()
        if isinstance(self.reff, shdom.GridDataEstimator):
            derivatives['reff'] = self.init_re_derivative()
        if isinstance(self.veff, shdom.GridPhaseEstimator):
            derivatives['veff'] = self.init_ve_derivative()
        return derivatives        
    
    def init_lwc_derivative(self):
        """TODO"""
        derivative = OrderedDict()
        derivative['lwc'] = shdom.GridData(self.lwc.grid, np.ones_like(self.lwc.data))
        derivative['albedo'] = shdom.GridData(self.grid, np.zeros(self.grid.shape))
        return derivative
        
        
    def init_re_derivative(self):
        """TODO"""
        derivatives = OrderedDict()
        for wavelength, mie in self.mie.items():
            extinct = mie.extinct.reshape((mie.size_distribution.nretab, mie.size_distribution.nvetab), order='F')
            ssalb = mie.ssalb.reshape((mie.size_distribution.nretab, mie.size_distribution.nvetab), order='F')
         
            dre = np.diff(mie.size_distribution.reff)        
            dextinct = np.diff(extinct, axis=0) / dre[:,np.newaxis]
            dssalb = np.diff(ssalb, axis=0) / dre[:,np.newaxis]
            dlegcoef = np.diff(mie.legcoef_2d, axis=-2) / dre[:,np.newaxis]
        
            # Define a derivative Mie object, last derivative is duplicated
            derivative = copy.deepcopy(mie)
            derivative._extinct = np.vstack((dextinct, dextinct[-1])).ravel(order='F')
            derivative._ssalb = np.vstack((dssalb, dssalb[-1])).ravel(order='F')
            if mie.table_type == 'SCALAR':
                derivative.legcoef = np.concatenate((dlegcoef, dlegcoef[:,-1][:,None]), axis=-2).reshape((mie.maxleg+1, -1), order='F')
            elif mie.table_type == 'VECTOR':
                derivative.legcoef = np.concatenate((dlegcoef, dlegcoef[:,-1][:,None]), axis=-2).reshape((6, mie.maxleg+1, -1), order='F')
    
            derivative.init_intepolators()
            derivatives[float_round(wavelength)] = derivative
        return derivatives         
       
   
    def init_ve_derivative(self):
        """TODO"""
        derivatives = OrderedDict()
        for wavelength, mie in self.mie.items():        
            extinct = self.extinct.reshape((self.size_distribution.nretab, self.size_distribution.nvetab), order='F')
            ssalb = self.ssalb.reshape((self.size_distribution.nretab, self.size_distribution.nvetab), order='F')
            if self.table_type == 'SCALAR':
                legcoef = self.legcoef.reshape((self.maxleg+1, self.size_distribution.nretab, self.size_distribution.nvetab), order='F')
            elif self.table_type == 'VECTOR':
                legcoef = self.legcoef.reshape((-1, self.maxleg+1, self.size_distribution.nretab, self.size_distribution.nvetab), order='F')
        
            dve = np.diff(self.size_distribution.veff)        
            dextinct = np.diff(extinct, axis=1) / dve[:,np.newaxis]
            dssalb = np.diff(ssalb, axis=1) / dve[:,np.newaxis]
            dlegcoef = np.diff(legcoef, axis=-1) / dve[:,np.newaxis]
        
            # Define a derivative Mie object, last derivative is duplicated
            derivative = copy.deepcopy(self)
            derivative._extinct = np.hstack((dextinct, dextinct[-1])).ravel(order='F')
            derivative._ssalb = np.hstack((dssalb, dssalb[-1])).ravel(order='F')
            if self.table_type == 'SCALAR':
                derivative.legcoef = np.concatenate((dlegcoef, dlegcoef[:,-1][:,None]), axis=-1).reshape((self.maxleg+1, -1), order='F')
            elif self.table_type == 'VECTOR':
                derivative.legcoef = np.concatenate((dlegcoef, dlegcoef[:,-1][:,None]), axis=-1).reshape((6, self.maxleg+1, -1), order='F')
        
            derivative.init_intepolators()
            derivatives[float_round(wavelength)] = derivative
        return derivatives  

       
    def get_lwc_derivative(self, wavelength):
        """TODO"""
        index  = shdom.GridData(self.grid, np.ones(self.grid.shape, dtype=np.int32))
        legen_table = shdom.LegendreTable(np.zeros((self.mie[float_round(wavelength)].maxleg+1), dtype=np.float32), 
                                          table_type=self.mie[float_round(wavelength)].table_type)       
        derivative = self.derivatives['lwc']
        scatterer = shdom.OpticalScattererDerivative(
            wavelength, 
            extinction=self.mie[float_round(wavelength)].get_extinction(derivative['lwc'], self.reff, self.veff),
            albedo=derivative['albedo'],
            phase=shdom.GridPhase(legen_table, index)) 
        return scatterer
    
    def get_phase_derivative(self, derivative, wavelength):
        """TODO"""
        scatterer = shdom.OpticalScattererDerivative(
            wavelength, 
            extinction=derivative[float_round(wavelength)].get_extinction(self.lwc, self.reff, self.veff),
            albedo=derivative[float_round(wavelength)].get_albedo(self.reff, self.veff),
            phase=derivative[float_round(wavelength)].get_phase(self.reff, self.veff)) 
        return scatterer  
    
    

    def get_derivative(self, derivative_type, wavelength):
        """TODO"""
        if derivative_type == 'lwc':
            derivative = self.get_lwc_derivative(wavelength)
        elif derivative_type == 'reff' or derivative_type == 'veff':
            derivative = self.get_phase_derivative(self.derivatives[derivative_type], wavelength)
        else:
            raise AttributeError('derivative type {} not supported'.format(derivative_type))
        return derivative
        

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
        for (name, estimator), state in zip(self.estimators.items(), states):
            estimator.set_state(state)
            self.scatterers[name] = estimator
    
    
    def get_state(self):
        """TODO"""
        state = np.empty(shape=(0),dtype=np.float64)
        for estimator in self.estimators.values():
            state = np.concatenate((state, estimator.get_state()))
        return state


    def get_bounds(self):
        """TODO
        bounds to a list of bounds that the scipy minimize expects."""
        bounds = []
        for estimator in self.estimators.values():
            bounds.extend(estimator.get_bounds())
        return bounds
    
    
    def get_derivatives(self, rte_solver):
        """TODO"""
        dext = np.zeros(shape=[rte_solver._nbpts, self.num_derivatives], dtype=np.float32)
        dalb = np.zeros(shape=[rte_solver._nbpts, self.num_derivatives], dtype=np.float32)
        diphase = np.zeros(shape=[rte_solver._nbpts, self.num_derivatives], dtype=np.float32)
    
        i=0
        for estimator in self.estimators.values():
            for dtype in estimator.derivatives.keys():
                derivative = estimator.get_derivative(dtype, rte_solver.wavelength)       
                resampled_derivative = derivative.resample(self.grid)
                extinction = resampled_derivative.extinction.data
                albedo = resampled_derivative.albedo.data
                iphase = resampled_derivative.phase.iphasep                
                
                # TODO: Check this
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
            negcheck=False,
            nscatangle=rte_solver._nscatangle,
            numphase=dnumphase,
            ml=rte_solver._ml,
            nlm=rte_solver._nlm,
            nleg=rte_solver._nleg,
            legen=dleg,
            deltam=False
        )                
        return dext, dalb, diphase, dleg, dphasetab, dnumphase
              
                 
    def compute_direct_derivative(self, rte_solver):
        """TODO"""
        
        if isinstance(rte_solver, shdom.RteSolverArray):
            rte_solver = rte_solver[0]
            
        self._direct_derivative_path, self._direct_derivative_ptr = \
            core.make_direct_derivative(
                dirflux=rte_solver._dirflux, 
                extdirp=rte_solver._pa.extdirp,                
                npts=rte_solver._npts,
                bcflag=rte_solver._bcflag,
                ipflag=rte_solver._ipflag,
                deltam=rte_solver._deltam,
                ml=rte_solver._ml,
                nleg=rte_solver._nleg,
                solarflux=rte_solver._solarflux,
                solarmu=rte_solver._solarmu,
                solaraz=rte_solver._solaraz,
                gridpos=rte_solver._gridpos,
                npx=rte_solver._pa.npx,
                npy=rte_solver._pa.npy,
                npz=rte_solver._pa.npz,
                numphase=rte_solver._pa.numphase,
                delx=rte_solver._pa.delx,
                dely=rte_solver._pa.dely,
                xstart=rte_solver._pa.xstart,
                ystart=rte_solver._pa.ystart,
                zlevels=rte_solver._pa.zlevels,
                tempp=rte_solver._pa.tempp,
                extinctp=rte_solver._pa.extinctp,
                albedop=rte_solver._pa.albedop,
                legenp=rte_solver._pa.legenp,
                iphasep=rte_solver._pa.iphasep,
                nzckd=rte_solver._pa.nzckd,
                zckd=rte_solver._pa.zckd,
                gasabs=rte_solver._pa.gasabs
            )       
        
    def core_grad(self, rte_solver, projection, radiance):
        """
        TODO
        """
        if isinstance(projection.npix, list):
            total_pix = np.sum(projection.npix)
        else:
            total_pix = projection.npix

        gradient, loss, radiance = core.gradient(
            dpath=self._direct_derivative_path, 
            dptr=self._direct_derivative_ptr,
            npx=rte_solver._pa.npx,
            npy=rte_solver._pa.npy,
            npz=rte_solver._pa.npz,
            delx=rte_solver._pa.delx,
            dely=rte_solver._pa.dely,                
            xstart=rte_solver._pa.xstart,
            ystart=rte_solver._pa.ystart,
            zlevels=rte_solver._pa.zlevels, 
            extdirp=rte_solver._pa.extdirp,
            uniformzlev=rte_solver._uniformzlev,
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
            solarflux=rte_solver._solarflux,
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
                negcheck=True,
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
        loss = np.sum(list(map(lambda x: x[1], output)))
        
        # Project rte_solver base grid gradient to the state space
        gradient = sum(list(map(lambda x: x[0], output)))
        gradient = gradient.reshape(self.grid.shape + tuple([self.num_derivatives]))
        images = sensor.make_images(np.concatenate(list(map(lambda x: x[2], output))), projection, num_channels)
        
        state_gradient = np.empty(shape=(0), dtype=np.float64)
        for i, estimator in enumerate(self.estimators.values()):
            state_gradient = np.concatenate(
                (state_gradient, estimator.project_gradient(GridData(self.grid, gradient[...,i])))
            )
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
        self._ckpt_periods = []
        self._ckpt_times = []
        self._callback_fns = []
        self._optimizer = None
        
        
    def attach_optimizer(self, optimizer):
        """TODO"""
        self._optimizer = optimizer
    
        
    def monitor_loss(self, ckpt_period=-1):
        """
        Monitor the loss.
        
        Parameters
        ----------
        ckpt_period: float
           time [seconds] between updates. setting ckpt_period=-1 will log at every iteration.
        """
        self._ckpt_periods.append(ckpt_period)
        self._ckpt_times.append(time.time())
        self._callback_fns.append(self.loss_cbfn) 
        

    def save_checkpoints(self, ckpt_period=-1):
        """Save a checkpoint of the Optimizer
        
        Parameters
        ----------
        ckpt_period: float
           time [seconds] between updates. setting ckpt_period=-1 will log at every iteration.
        """
        self._ckpt_periods.append(ckpt_period)
        self._ckpt_times.append(time.time())        
        self._callback_fns.append(self.save_ckpt_cbfn)
            
        
    def monitor_shdom_iterations(self, ckpt_period=-1):
        """Monitor the number of SHDOM iterations.
        
        Parameters
        ----------
        ckpt_period: float
           time [seconds] between updates. setting ckpt_period=-1 will log at every iteration.
        """
        self._ckpt_periods.append(ckpt_period)
        self._ckpt_times.append(time.time())
        self._callback_fns.append(self.shdom_iterations_cbfn)         
        
        
    def monitor_scatterer_error(self, estimator_name, ground_truth, ckpt_period=-1):
        """
        Monitor relative and overall mass error (epsilon, delta) as defined at:
          Amit Aides et al, "Multi sky-view 3D aerosol distribution recovery".
        
        Parameters
        ----------
        estimator_name: str
            The name of the scatterer to monitor
        ground_truth: shdom.Scatterer
            The ground truth medium.
        ckpt_period: float
           time [seconds] between updates. setting ckpt_period=-1 will log at every iteration.
        """
        self._ckpt_periods.append(ckpt_period)
        self._ckpt_times.append(time.time())
        self._callback_fns.append(self.scatterer_error_cbfn)
        if hasattr(self, '_ground_truth'):
            self._ground_truth[estimator_name] = ground_truth
        else:
            self._ground_truth = OrderedDict({estimator_name: ground_truth})
            
    def monitor_domain_mean(self, estimator_name, ground_truth, ckpt_period=-1):
        """
        Monitor domain mean and compare to ground truth over iterations.

        Parameters
        ----------
        estimator_name: str
            The name of the scatterer to monitor
        ground_truth: shdom.Scatterer
            The ground truth medium.
        ckpt_period: float
           time [seconds] between updates. setting ckpt_period=-1 will log at every iteration.
        """
        self._ckpt_periods.append(ckpt_period)
        self._ckpt_times.append(time.time())
        self._callback_fns.append(self.domain_mean_cbfn)
        if hasattr(self, '_ground_truth'):
            self._ground_truth[estimator_name] = ground_truth
        else:
            self._ground_truth = OrderedDict({estimator_name: ground_truth})
        
        
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
        num_images = len(acquired_images)
        self._ckpt_periods.append(ckpt_period)
        self._ckpt_times.append(time.time())
        self._callback_fns.append(self.estimated_images_cbfn)
        self._image_vmax = [image.max() * 1.25 for image in acquired_images]
        self._image_titles = ['Retrieval Image {}'.format(view) for view in range(num_images)]
        acq_titles = ['Acquiered Image {}'.format(view) for view in range(num_images)]
        self.write_image_list(0, acquired_images, acq_titles, vmax=self._image_vmax)


    def save_ckpt_cbfn(self):
        timestr = time.strftime("%H%M%S")
        path = os.path.join(self.writer._tf_writer.log_dir,  timestr + '.ckpt')
        self.optimizer.save(path)
        
    def loss_cbfn(self):
        """Callback function that is called (every optimizer iteration) for loss monitoring."""
        self.tf_writer.add_scalar('loss', self.optimizer.loss, self.optimizer.iteration)    
    
    def estimated_images_cbfn(self):
        """Callback function the is called every optimizer iteration image monitoring is set."""
        self.write_image_list(self.optimizer.iteration, self.optimizer.images, self._image_titles, self._image_vmax)
    
    def shdom_iterations_cbfn(self):
        """Callback function that is called (every optimizer iteration) for shdom iteration monitoring"""
        self.tf_writer.add_scalar('total shdom iterations', self.optimizer.rte_solver.num_iterations, self.optimizer.iteration)       
        
    def scatterer_error_cbfn(self):
        """Callback function for monitoring parameter error measures."""
        for scatterer_name, gt_scatterer in self._ground_truth.items():
            est_scatterer = self.optimizer.medium.get_scatterer(scatterer_name)
            for parameter_name, parameter in est_scatterer.estimators.items():
                est_param = parameter.data[parameter.mask.data].ravel()
                gt_param = getattr(gt_scatterer, parameter_name).data[parameter.mask.data].ravel()
                delta = (np.linalg.norm(est_param, 1) - np.linalg.norm(gt_param, 1)) / np.linalg.norm(gt_param, 1)
                epsilon = np.linalg.norm((est_param - gt_param), 1) / np.linalg.norm(gt_param,1)
                self.tf_writer.add_scalar(scatterer_name + ' delta ' + parameter_name, delta, self.optimizer.iteration)
                self.tf_writer.add_scalar(scatterer_name + ' epsilon ' + parameter_name, epsilon, self.optimizer.iteration)           
        
    def domain_mean_cbfn(self):
        """Callback function for monitoring parameter error measures."""
        for scatterer_name, gt_scatterer in self._ground_truth.items():
            est_scatterer = self.optimizer.medium.get_scatterer(scatterer_name)
            for parameter_name, parameter in est_scatterer.estimators.items():
                est_param = parameter.data[parameter.mask.data].mean()
                gt_param = getattr(gt_scatterer, parameter_name).data[parameter.mask.data].mean()
                self.tf_writer.add_scalars(
                    main_tag=scatterer_name + ' mean ' + parameter_name,
                    tag_scalar_dict={'esimated': est_param, 'true': gt_param}, 
                    global_step=self.optimizer.iteration
                )
                
     
    def write_image_list(self, global_step, images, titles, vmax=None):
        """
        Write an image list to tensorboardX.
    
        Parameters
        ----------
        global_step: integer,
            The global step of the optimizer.
        images: list
            List of images to be logged onto tensorboard.
        titles: list
            List of strings that will title the corresponding images on tensorboard.
        vmax: list or scalar, optional
            List or a single of scaling factor for the image contrast equalization
        """
    
        if np.isscalar(vmax) or vmax is None:
            vmax = [vmax]*len(images)        
    
        assert len(images) == len(titles), 'len(images) != len(titles): {} != {}'.format(len(images), len(titles))
        assert len(vmax) == len(titles), 'len(vmax) != len(images): {} != {}'.format(len(vmax), len(times))
    
        for image, title, vm in zip(images, titles, vmax):
            # for polychromatic
            if image.ndim == 3:
                self.tf_writer.add_images(
                    tag=title, 
                    img_tensor=(np.repeat(np.expand_dims(image, 2), 3, axis=2) / vm),
                    dataformats='HWCN',
                    global_step=global_step
                )
            # for monochromatic
            else:
                self.tf_writer.add_image(
                    tag=title, 
                    img_tensor=(image / vm),
                    dataformats='HW',
                    global_step=global_step
                ) 
                
    def check_update_time(self, index):
        """
        Check if it is time to update the callback function.
        
        Parameters
        ----------
        index: integer,
            The index of the callback function in the list.
        """
        time_passed = time.time() - self.ckpt_times[index] 
        if time_passed > self.ckpt_periods[index]:
            self._ckpt_times[index] = time.time()
            return True
        else:
            return False
        
        
    @property
    def callback_fns(self):
        return self._callback_fns
    
    @property
    def ckpt_periods(self):
        return self._ckpt_periods
    
    @property
    def ckpt_times(self):
        return self._ckpt_times
    
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
        A threshold on radiances is used to produce a mask and preform space carving.
        
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
        gradient, loss, images = self.medium.compute_gradient(
            rte_solver=self.rte_solver,
            measurements=self.measurements,
            n_jobs=self._n_jobs
        )
        self._loss = loss
        self._images = images
        return loss, gradient
    

            
    def callback(self, state):
        """
        The callback function invokes the callbacks defined by the writer (if any). 
        Additionally it keeps track of the iteration number 
        and if the checkpoint period has passed it saves a checkpoint of the Optimizer
        """
        self._iteration += 1
        
        # Writer callback functions
        if self.writer is not None:
            for index, callbackfn in enumerate(self.writer.callback_fns):
                if self.writer.check_update_time(index) == True:
                    callbackfn()
        

    def minimize(self, options, method='L-BFGS-B', n_jobs=1):
        """
        Minimize the cost function with respect to the parameters defined.
        
        Parameters
        ----------
        options: Dict
            The option dictionary
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
        self.rte_solver.set_medium(self.medium)
        self.rte_solver.init_solution()
        self.medium.compute_direct_derivative(self.rte_solver)
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
        self.rte_solver.set_medium(self.medium)
        self.rte_solver.make_direct()
        self.rte_solver.solve(maxiter=100, verbose=False)        
        
        
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
        params['rte_param_dict'] = rte_solver.get_param_dict()
        file = open(path,'wb')
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
        file = open(path,'rb')
        data = file.read()
        file.close()        
        params = pickle.loads(data)
        
        # Create an RTE solver object
        rte_param_dict = params.pop('rte_param_dict')
        if len(rte_param_dict['solver_parameters']) == 1:
            rte_solver = shdom.RteSolver()
        else:
            rte_solver = shdom.RteSolverArray()
        rte_solver.set_param_dict(rte_param_dict)         
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