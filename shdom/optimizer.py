"""
TODO
"""

import shdom
import numpy as np
import core 
from scipy.optimize import minimize
from shdom import GridData


class Optimizer(object):
    """
    TODO
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
        
    def set_measurements(self, measurements):
        """TODO"""
        self._measurements = measurements
        
    def set_rte_solver(self, rte_solver):
        """TODO"""
        self._rte_solver = rte_solver

    def set_cloud_mask(self, cloud_mask):
        """TODO"""
        self._cloud_mask = cloud_mask  
    
    def set_known_medium(self, medium):
        """TODO"""
        self._known_medium = medium

    def add_parameter(self, parameter):
        """TODO"""
        if parameter.__class__ is not shdom.parameters.Extinction:
            raise NotImplementedError('Only Extinction optimization is supported.')
        
        self._parameters.append(parameter)
        self._extinction_dependency = self.extinction_dependency or parameter.extinction_dependency
        self._albedo_dependency = self.albedo_dependency or parameter.albedo_dependency
        self._phase_dependency = self.phase_dependency or parameter.phase_dependency
        
    def update_rte_solver(self, state):
        """
        TODO
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
        Gradient and cost computation by ray tracing the domain (see shdomsub4.f)
        
        
        Returns
        -------
        
        Notes
        -----
        """ 
        gradient, cost = core.ext_gradient(
            nx=self.rte_solver._nx,
            ny=self.rte_solver._ny,
            nz=self.rte_solver._nz,
            bcflag=self.rte_solver._bcflag,
            ipflag=self.rte_solver._ipflag,   
            npts=self.rte_solver._npts,
            nbpts=self._rte_solver._nbpts,
            ncells=self.rte_solver._ncells,
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
            npix=self.sensors.npix,          
            srctype=self.rte_solver._srctype,
            sfctype=self.rte_solver._sfctype,
            units=self.rte_solver._units,
            measurements=self.radiances,
            rshptr=self.rte_solver._rshptr,
            radiance=self.rte_solver._radiance
        )  
        gradout = self.basegrid_projection(gradient)       
        return gradout, cost
    
    
    def objective_fun(self, state):
        """TODO"""
        self.update_rte_solver(state)     
        self.rte_solver.solve(maxiter=100, verbose=False)  
        return self.extinction_gradient_cost()[1]
        
        
    def gradient(self, state):
        """TODO"""
        self.update_rte_solver(state)
        parameter_gradient = [self.extinction_gradient_cost()[0]]
        state_gradient = np.concatenate(map(lambda grad: grad[self.gradient_mask].ravel(), parameter_gradient))
        return state_gradient
    
    
    def minimize(self, options, method='L-BFGS-B'):
        """
        TODO
        """
        
        if method != 'L-BFGS-B':
            raise NotImplementedError('Optimization method not implemented')
        
        # Init masks for parameters
        self._num_parameters = []
        for param in self.parameters:
            if self.cloud_mask:
                param.set_mask(self.cloud_mask)
                self._num_parameters.append(param.num_parameters)
                
        # Init mask for gradient 
        if self.cloud_mask:
            self._gradient_mask = self.cloud_mask.data.ravel()
            if self.known_medium:
                grid = self.known_medium.grid + self.cloud_mask.grid
                self._gradient_mask = np.array(self.cloud_mask.resample(grid, method='nearest'), dtype=np.bool).ravel()
            
        initial_state = self.parameters_to_state()
        result = minimize(fun=self.objective_fun, 
                          x0=initial_state, 
                          method=method, 
                          jac=self.gradient,
                          bounds=self.get_bounds(),
                          options=options)
        return result
    
    
    def parameters_to_state(self):
        """TODO"""
        state = np.concatenate(map(lambda param: param.data[param.mask].ravel(), self.parameters))
        return state
    
    
    def get_bounds(self):
        """TODO"""
        bounds = map(lambda param: param.bounds, self.parameters)
        return bounds


    def set_state(self, state):
        """TODO"""
        state = np.concatenate(map(lambda param: param.data, self.parameters))
        return state    


    def basegrid_projection(self, field):
        """TODO: Check against fortran subroutine"""
        for icell in range(self.rte_solver._nbcells, self.rte_solver._ncells):
            # Find base cell for icell
            base_cell = icell
            while (self.rte_solver._treeptr[0, base_cell] > 0):
                base_cell = self.rte_solver._treeptr[0, base_cell] - 1
    
            # The two most distant points of the 8 define the interpolation
            base_point0 = self.rte_solver._gridptr[0, base_cell] - 1
            base_point7 = self.rte_solver._gridptr[7, base_cell] - 1
            bx0, by0, bz0 = self.rte_solver._gridpos[:, base_point0]
            bx7, by7, bz7 = self.rte_solver._gridpos[:, base_point7]
            dx, dy, dz = bx7-bx0, by7-by0, bz7-bz0
    
            # loop over gridpoints belonging to icell and trilin interpolate
            for gridpoint in self.rte_solver._gridptr[:, icell]-1:
                if (field[gridpoint] != 0) & (gridpoint >= self.rte_solver._nbpts):
                    u, v, w = 0, 0, 0
                    if dx>0.0:
                        u = (bx7 - self.rte_solver._gridpos[0, gridpoint])/dx
                    if dy>0.0:
                        v = (by7 - self.rte_solver._gridpos[1, gridpoint])/dy
                    if dz>0.0:
                        w = (bz7 - self.rte_solver._gridpos[2, gridpoint])/dz
    
                    f = np.empty(8)
                    f[7] = (1-u)*(1-v)*(1-w)
                    f[6] =   u  *(1-v)*(1-w)
                    f[5] = (1-u)*  v  *(1-w)
                    f[4] =   u  *  v  *(1-w)
                    f[3] = (1-u)*(1-v)*  w
                    f[2] =   u  *(1-v)*  w
                    f[1] = (1-u)*  v  *  w
                    f[0] =   u  *  v  *  w
    
                    for bgridpoint, _f in zip(self.rte_solver._gridptr[:,base_cell]-1, f):
                        field[bgridpoint] += _f*field[gridpoint]
                        field[gridpoint]   = 0
        return field[:self.rte_solver._nbpts] 
        
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