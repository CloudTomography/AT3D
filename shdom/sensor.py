"""
TODO: description
"""

import core
import numpy as np
from shdom import BoundingBox
import itertools

norm = lambda x: x / np.linalg.norm(x, axis=0)

class Sensor(object):
    """
    TODO: add description.
    Abstract Sensor class to be inherited by the different types of sensors.
    Each sensor internally defines an arrays of pixel locations (x,y,z) in km and directions (phi, mu).
    """
    def __init__(self):
        self._x = None
        self._y = None
        self._z = None
        self._mu = None
        self._phi = None
        self._npix = None
        self._resolution = None
        self.type = 'AbstractSensor'
        
        
    def project(self, projection_matrix, point_array):
        """TODO"""
        homogenic_point_array = np.pad(point_array,((0,1),(0,0)),'constant', constant_values=1)
        return np.dot(projection_matrix, homogenic_point_array)    


    def render(self, rte_solver):
        """TODO"""
        visout = core.render(
            camx=self.x,
            camy=self.y,
            camz=self.z,
            cammu=self.mu,
            camphi=self.phi,
            npix=self._npix,             
            nx=rte_solver._nx,
            ny=rte_solver._ny,
            nz=rte_solver._nz,
            bcflag=rte_solver._bcflag,
            ipflag=rte_solver._ipflag,   
            npts=rte_solver._npts,
            ncells=rte_solver._ncells,
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
            iphase=rte_solver._iphase,
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
            extinct=rte_solver._extinct,
            albedo=rte_solver._albedo,
            legen=rte_solver._legen,            
            dirflux=rte_solver._dirflux,
            fluxes=rte_solver._fluxes,
            source=rte_solver._source,          
            srctype=rte_solver._srctype,
            sfctype=rte_solver._sfctype,
            units=rte_solver._units
        )
        return visout.reshape(self.resolution, order='F')      

    @property
    def x(self):
        return self._x
    
    @property
    def y(self):
        return self._y
    
    @property
    def z(self):
        return self._z
    
    @property
    def mu(self):
        return self._mu
    
    @property
    def phi(self):
        return self._phi
    
    @property
    def zenith(self):
        return np.rad2deg(np.arccos(mu))
 
    @property
    def azimuth(self):
        return np.rad2deg(self.phi)
    
    @property 
    def resolution(self):
        return self._resolution
    
    @property 
    def npix(self):
        return self._npix    
    
    
class OrthographicSensor(Sensor):
    """
    TODO: parallel rays description 
    """
    
    def __init__(self, bounding_box, x_resolution, y_resolution, azimuth, zenith, altitude='TOA'):
        super(OrthographicSensor, self).__init__()
        self.type = 'OrthographicSensor'
        
        self._x_resolution = x_resolution
        self._y_resolution = y_resolution
        
        mu = np.cos(np.deg2rad(zenith))
        phi = np.deg2rad(azimuth)
        if altitude == 'TOA':
            self._altitude = bounding_box.zmax
        else:
            assert (type(altitude) == float or type(altitude) == int), 'altitude of incorrect type'
            self._altitude = altitude        
            
        # Project the bounding box onto the image plane
        alpha = np.sqrt(1 - mu**2) * np.cos(phi) / mu
        beta  = np.sqrt(1 - mu**2) * np.sin(phi) / mu
        projection_matrix = np.array([
            [1, 0, -alpha, alpha* self.altitude],
            [0, 1,  -beta, beta * self.altitude],
            [0, 0,      0,        self.altitude]
        ])
        
        self.projection_matrix = projection_matrix
        bounding_box_8point_array = np.array(list(itertools.product([bounding_box.xmin, bounding_box.xmax], 
                                                                    [bounding_box.ymin, bounding_box.ymax],
                                                                    [bounding_box.zmin, bounding_box.zmax]))).T
        projected_bounding_box = self.project(projection_matrix, bounding_box_8point_array)
            
        # Use projected bounding box to define image sampling
        x_s, y_s = projected_bounding_box[:2,:].min(axis=1)
        x_e, y_e = projected_bounding_box[:2,:].max(axis=1)
        x = np.arange(x_s, x_e+1e-6, self.x_resolution)
        y = np.arange(y_s, y_e+1e-6, self.y_resolution)
        z = self.altitude
        self._x, self._y, self._z, self._mu, self._phi = np.meshgrid(x, y, z, mu, phi)
        self._x = self._x.ravel().astype(np.float32)
        self._y = self._y.ravel().astype(np.float32)
        self._z = self._z.ravel().astype(np.float32)
        self._mu = self._mu.ravel().astype(np.float32)
        self._phi = self._phi.ravel().astype(np.float32)
        self._npix = self.x.size
        self._resolution = [x.size, y.size]

    
    @property 
    def altitude(self):
        return self._altitude
    
    @property 
    def x_resolution(self):
        return self._x_resolution
     
    @property 
    def y_resolution(self):
        return self._y_resolution
    
    

class ProjectiveSensor(Sensor):
    """
    TODO: description 
    
    notes:
    -----
    pinhole model
    """
    
    def __init__(self, fov, nx, ny, x, y, z):
        super(ProjectiveSensor, self).__init__()
        self.type = 'ProjectiveSensor'

        self._resolution = [nx, ny]
        self._npix = nx*ny
        self._position = np.array([x, y, z], dtype=np.float32)
        self._fov = fov
        self._focal = 1.0 / np.tan(np.deg2rad(fov) / 2.0)
        self._k = np.array([[self._focal, 0, 0],
                            [0, self._focal, 0],
                            [0, 0, 1]], dtype=np.float32)
        self._inv_k = np.linalg.inv(self._k)
        self._rotation_matrix = np.eye(3)
        x_c, y_c, z_c = np.meshgrid(np.linspace(-1, 1, nx), np.linspace(-1, 1, ny), 1.0)        
        self._homogeneous_coordinates = np.stack([x_c.ravel(), y_c.ravel(), z_c.ravel()])
        self.update_global_coordinates()
        
        
    def update_global_coordinates(self):
        """TODO"""
        x_c, y_c, z_c = norm(np.matmul(
            self._rotation_matrix, np.matmul(self._inv_k, self._homogeneous_coordinates)))

        self._mu = -z_c.astype(np.float32)
        self._phi = (np.arctan2(y_c, x_c) + np.pi).astype(np.float32)
        self._x = np.full(self.npix, self.position[0], dtype=np.float32)
        self._y = np.full(self.npix, self.position[1], dtype=np.float32)
        self._z = np.full(self.npix, self.position[2], dtype=np.float32)
    
    
    def look_at_transform(self, point, up):
        """TODO"""
        up = np.array(up)
        direction = np.array(point) - self.position
        zaxis = norm(direction)
        xaxis = norm(np.cross(up, zaxis))
        yaxis = np.cross(zaxis, xaxis)
        self._rotation_matrix = np.stack((xaxis, yaxis, zaxis), axis=1)
        self.update_global_coordinates()
        
        
    def rotate_transform(self, axis, angle):
        """TODO"""
        
        assert axis in ['x', 'y', 'z'], 'axis parameter can only recieve "x", "y" or "z"'
        angle = np.deg2rad(angle)
        if axis == 'x':
            rot = np.array([[1, 0, 0],
                            [0, np.cos(angle), -np.sin(angle)],
                            [0, np.sin(angle), np.cos(angle)]], dtype=np.float32)
        elif axis == 'y':
            rot = np.array([[np.cos(angle), 0, np.sin(angle)],
                            [0, 1, 0],
                            [-np.sin(angle), 0, np.cos(angle)]], dtype=np.float32)
        elif axis == 'z':
            rot = np.array([[np.cos(angle), -np.sin(angle), 0],
                            [np.sin(angle), np.cos(angle), 0],
                            [0, 0, 1]], dtype=np.float32)        
        
        self._rotation_matrix = np.matmul(self._rotation_matrix, rot)
        self.update_global_coordinates()
        
        
    def plot(self, ax, xlim, ylim, zlim, length=0.1):
        mu = -self.mu.reshape(self.resolution)[[0, -1, 0, -1],[0, 0, -1, -1]]
        phi = np.pi + self.phi.reshape(self.resolution)[[0, -1, 0, -1],[0, 0, -1, -1]]
        u = np.sqrt(1 - mu**2) * np.cos(phi)
        v = np.sqrt(1 - mu**2) * np.sin(phi)
        w = mu
        x = np.full(4, self.position[0], dtype=np.float32)
        y = np.full(4, self.position[1], dtype=np.float32)
        z = np.full(4, self.position[2], dtype=np.float32)
        ax.set_aspect('equal')
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_zlim(*zlim)
        ax.quiver(x, y, z, u, v, w, length=length, pivot='tail')
        
        
    @property
    def position(self):
        return self._position
