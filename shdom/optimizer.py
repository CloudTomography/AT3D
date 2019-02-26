import shdom
import numpy as np
import core 
from scipy.optimize import minimize


class Optimizer(object):
    """
    TODO
    """
    def __init__(self):
        pass
    
    
    def set_measurements(self, measurements):
        self._measurements = np.array(measurements, dtype=np.float32).ravel()
        
    def set_rte_solver(self, rte_solver):
        self._rte_solver = rte_solver
        
        
    def set_sensor(self, sensor):
        self._sensor = sensor


  
    def compute_gradient_cost(self, x):
        """
        Gradient and cost computation by ray tracing the domain (see shdomsub4.f)
        
        Parameters
        ----------
        
        Returns
        -------
        
        Notes
        -----
        """ 
        ext = np.zeros((3,3,3), dtype=np.float32)
        ext[1,1,1] = x
        self.rte_solver._pa.extinctp = ext.ravel()
        gradient, cost = core.compute_gradient(
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
            camx=self.sensor.x,
            camy=self.sensor.y,
            camz=self.sensor.z,
            cammu=self.sensor.mu,
            camphi=self.sensor.phi,
            npix=self.sensor.npix,          
            srctype=self.rte_solver._srctype,
            sfctype=self.rte_solver._sfctype,
            units=self.rte_solver._units,
            measurements=self.measurements,
            rshptr=self.rte_solver._rshptr,
            radiance=self.rte_solver._radiance
        )  
        gradout = np.array(self.basegrid_projection(gradient).reshape(3,3,3)[1,1,1])
        return gradout, cost
    
    
    def cost(self, x):
        ext = np.zeros((3,3,3), dtype=np.float32)
        ext[1,1,1] = x
        self.rte_solver._pa.extinctp = ext.ravel()        
        self.rte_solver.solve(maxiter=100, verbose=False)  
        return self.compute_gradient_cost(x)[1]
        
        
    def minimize(self, init, options, method='L-BFGS-B'):
        """
        TODO
        """
        
        # Initialize the rte_solver to the initial medium
        self.rte_solver.init_medium(init)
        
        if method != 'L-BFGS-B':
            raise NotImplementedError('Optimization method not implemented')

        gradient = lambda x: self.compute_gradient_cost(x)[0]
        
        x0 = np.array(0.01, dtype=np.float64)
        bounds = [(0.0, None)]
        
        result = minimize(fun=self.cost, 
                          x0=x0, 
                          method=method, 
                          jac=gradient,
                          bounds=bounds,
                          options=options)
        return result
    
    
    def basegrid_projection(self, field):
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
    def measurements(self):
        return self._measurements
    
    @property
    def sensor(self):
        return self._sensor