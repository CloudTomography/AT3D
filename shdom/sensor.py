"""
Camera, Sensor and Projection related objects used for rendering.
"""
import core
import numpy as np
import shdom 
import itertools
import dill as pickle
from joblib import Parallel, delayed
import copy

norm = lambda x: x / np.linalg.norm(x, axis=0)

    
class Sensor(object):
    """
    A sensor class to be inhireted by specific sensor types.
    A Sensor needs to define a render method.
    """
    def __init__(self):
        self._type = 'Sensor'
    
    def render(self, rte_solver, projection):
        """TODO"""
        
        if isinstance(projection.npix, list):
            total_pix = np.sum(projection.npix)
        else:
            total_pix = projection.npix 
            
        output = core.render(
            ncs=rte_solver._ncs,
            nstokes=rte_solver._nstokes,
            nstleg=rte_solver._nstleg,
            camx=projection.x,
            camy=projection.y,
            camz=projection.z,
            cammu=projection.mu,
            camphi=projection.phi,
            npix=total_pix,             
            nx=rte_solver._nx,
            ny=rte_solver._ny,
            nz=rte_solver._nz,
            bcflag=rte_solver._bcflag,
            ipflag=rte_solver._ipflag,   
            npts=rte_solver._npts,
            ncells=rte_solver._ncells,
            ml=rte_solver._ml,
            mm=rte_solver._mm,
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
            legen=rte_solver._legen,            
            dirflux=rte_solver._dirflux,
            fluxes=rte_solver._fluxes,
            source=rte_solver._source,          
            srctype=rte_solver._srctype,
            sfctype=rte_solver._sfctype,
            units=rte_solver._units,
            total_ext=rte_solver._total_ext[:rte_solver._npts],
            npart=rte_solver._npart)    
        
        return output
    
    
    @property
    def type(self):
        return self._type
    
    
class RadianceSensor(Sensor):
    """
    A Radiance sensor measures monochromatic radiances in [w/(m^2*sr*micron)]
    """
    def __init__(self):
        super(RadianceSensor, self).__init__()
        self._type = 'RadianceSensor'
        
    def render(self, rte_solver, projection, n_jobs=1, verbose=0):
        """
        The render method integrates a pre-computed in-scatter field (source function) J over the projection gemoetry.
        The source code for this function is in src/unoplarized/shdomsub4.f. 
        It is a modified version of the original SHDOM visualize_radiance subroutine in src/unpolarized/shdomsub2.f.
        
        If n_jobs>1 than parallel rendering is used where all pixels are distributed amongst all workers
        
        
        Parameters
        ----------
        rte_solver: shdom.RteSolver object
            The RteSolver with the precomputed radiative transfer solution (RteSolver.solve method).
        projection: shdom.Projection object
            The Projection specifying the sensor camera geomerty.
        n_jobs: int, default=1
            The number of jobs to divide the rendering into.
        verbose: int, default=0
            How much verbosity in the parallel rendering proccess.
            
        Returns
        -------
        radiance: np.array(shape=(projection.resolution), dtype=np.float)
            The rendered radiances.
        
        Notes
        -----
        For a small amout of pixels parallel rendering is slower due to communication overhead.
        """
        
        # If rendering several atmospheres (e.g. multi-spectral rendering)
        if isinstance(rte_solver, shdom.RteSolverArray):
            num_channels = rte_solver.num_solvers  
            rte_solvers = rte_solver
        else:
            num_channels = 1
            rte_solvers = [rte_solver]
            
        # Parallel rendering using multithreading (threadsafe Fortran)
        if n_jobs > 1:
            radiance = Parallel(n_jobs=n_jobs, backend="threading", verbose=verbose)(
                delayed(super(RadianceSensor, self).render, check_pickle=False)(
                    rte_solver=rte_solver,
                    projection=projection) for rte_solver, projection in 
                itertools.product(rte_solvers, projection.split(n_jobs)))  
            
        # Sequential rendering
        else:
            radiance = [super(RadianceSensor, self).render(rte_solver, projection) for rte_solver in rte_solvers]

        radiance = np.concatenate(radiance) 
        images = self.make_images(radiance, projection, num_channels)
        return images
            
        
    def make_images(self, radiance, projection, num_channels):
        """
        Split into Multiview, Multi-channel images (channel last)
        TODO
        """
        multiview = isinstance(projection, shdom.MultiViewProjection)
        multichannel = num_channels > 1
        
        if multichannel:
            radiance = np.array(np.split(radiance, num_channels)).T     
            
        if multiview: 
            split_indices = np.cumsum(projection.npix[:-1])        
            radiance = np.split(radiance, split_indices)
            
            if multichannel:
                radiance = [
                    image.reshape(resolution + [num_channels], order='F')
                    for image, resolution in zip(radiance, projection.resolution)
                ]
            else:
                radiance = [
                    image.reshape(resolution, order='F') 
                    for image, resolution in zip(radiance, projection.resolution) 
                ]                  
        else:
            new_shape = projection.resolution
            if multichannel:
                new_shape.append(num_channels)       
            radiance = radiance.reshape(new_shape, order='F') 
                
        return radiance         
    
 
class StokesSensor(Sensor):
    """
    A StokesSensor measures monochromatic stokes vector [I, U, Q, V].
    """
    def __init__(self):
        super(StokesSensor, self).__init__()
        self._type = 'StokesSensor'
    
    def render(self, rte_solver, projection, n_jobs=1, verbose=0):
        """      
        The render method integrates a pre-computed stokes vector in-scatter field (source function) J over the sensor gemoetry.
        The source code for this function is in src/polarized/shdomsub4.f. 
        It is a modified version of the original SHDOM visualize_radiance subroutine in src/polarized/shdomsub2.f.
        
        If n_jobs > 1 than parallel rendering is used where all pixels are distributed amongst all workers
        
        Parameters
        ----------
        rte_solver: shdom.RteSolver object
            The RteSolver with the precomputed radiative transfer solution (RteSolver.solve method).
        projection: shdom.Projection object
            The Projection specifying the sensor camera geomerty.
        n_jobs: int, default=1
            The number of jobs to divide the rendering into.
        verbose: int, default=0
            How much verbosity in the parallel rendering proccess.
            
            
        Returns
        -------
        stokes: np.array(shape=(nstokes, sensor.resolution), dtype=np.float32)
            The rendered radiances.

            
        Notes
        -----
        For a small amout of pixels parallel rendering is slower due to communication overhead.
        """
        multiview = isinstance(projection, shdom.MultiViewProjection)
        multichannel = isinstance(rte_solver, shdom.RteSolverArray)
        
        # If rendering several atmospheres (e.g. multi-spectral rendering)
        rte_solvers = rte_solver if multichannel else [rte_solver]
        
        # Parallel rendering using multithreading (threadsafe Fortran)
        if n_jobs > 1:
            stokes = Parallel(n_jobs=n_jobs, backend="threading", verbose=verbose)(
                delayed(super(StokesSensor, self).render, check_pickle=False)(
                    rte_solver=rte_solver,
                    projection=projection) for rte_solver, projection in 
                itertools.product(rte_solvers, projection.split(n_jobs)))  
            
        # Sequential rendering
        else:      
            stokes = [super(StokesSensor, self).render(rte_solver, projection) for rte_solver in rte_solvers]
          
        stokes = np.hstack(stokes) 
        images = make_images(stokes, projection, num_channels)
        return images
        
    def make_images(self, stokes, projection, num_channels):
        """
        TODO 
        Split into Multiview, Multi-channel images (channel last)
        """
        multiview = isinstance(projection, shdom.MultiViewProjection)
        multichannel = num_channels > 1        

        if multichannel:
            stokes = np.array(np.split(stokes, num_channels, axis=-1)).transpose([1, 2, 0]).squeeze()      
    
        if multiview: 
            split_indices = np.cumsum(projection.npix[:-1])        
            stokes = np.split(stokes, split_indices)
    
            if multichannel:
                stokes = [
                    image.reshape(resolution + [num_channels], order='F')
                    for image, resolution in zip(stokes, projection.resolution)
                ]
            else:
                stokes = [
                    image.reshape(resolution, order='F') 
                    for image, resolution in zip(stokes, projection.resolution) 
                ]                  
        else:
            new_shape =  [stokes.shape[0]] + projection.resolution
            if multichannel:
                new_shape.append(num_channels)       
            stokes = stokes.reshape(new_shape, order='F')        
            
        return stokes


class DolpAolpSensor(StokesSensor):
    """
    A DolpAolp measures monochromatic Degree and angle of Linear Polarization.
    """
    def __init__(self):
        super(StokesSensor, self).__init__()
        self._type = 'DolpAolpSensor'
    
    def render(self, rte_solver, projection, n_jobs=1, verbose=0):
        """   
        The render method integrates a pre-computed stokes vector in-scatter field (source function) J over the sensor gemoetry.
        The source code for this function is in src/polarized/shdomsub4.f. 
        It is a modified version of the original SHDOM visualize_radiance subroutine in src/polarized/shdomsub2.f.
        
        If n_jobs>1 than parallel rendering is used where all pixels are distributed amongst all workers
        
        Parameters
        ----------
        rte_solver: shdom.RteSolver object
            The RteSolver with the precomputed radiative transfer solution (RteSolver.solve method).
        projection: shdom.Projection object
            The Projection specifying the sensor camera geomerty.
        n_jobs: int, default=1
            The number of jobs to divide the rendering into.
        verbose: int, default=0
            How much verbosity in the parallel rendering proccess.
            
            
        Returns
        -------
        dolp: np.array(shape=(sensor.resolution), dtype=np.float32)
            Degree of Linear Polarization
        aolp: np.array(shape=(sensor.resolution), dtype=np.float32)
            Angle of Linear Polarization
        """
        stokes = super(DolpAolpSensor, self).render(rte_solver, projection, n_jobs, verbose)
        
        indices = stokes[0] > 0.0
        dolp = np.zeros_like(stokes[0])
        aolp = np.zeros_like(stokes[0])
        i, q, u = stokes[0][indices], stokes[1][indices], stokes[2][indices]
        dolp[indices] = np.sqrt(q**2 + u**2) / i
        aolp[indices] = (180.0/np.pi) * 0.5 * np.arctan2(u, q)
        
        # Choose the best range for the angle of linear polarization (-90 to 90 or 0 to 180)
        aolp1 = aolp.reshape(-1, aolp.shape[-1])
        aolp2 = aolp1.copy()
        aolp2[aolp2 < 0.0] += 180.0
        std1 = np.std(aolp1, axis=0)
        std2 = np.std(aolp2, axis=0)
        aolp[..., std2 < std1] = aolp2.reshape(aolp.shape)[..., std2 < std1]
        
        return dolp, aolp
    
    
class Projection(object):
    """
    Abstract Projection class to be inherited by the different types of projections.
    Each projection defines an arrays of pixel locations (x,y,z) in km and directions (phi, mu).
    """
    def __init__(self, x=None, y=None, z=None, mu=None, phi=None):
        self._x = x
        self._y = y
        self._z = z
        self._mu = mu
        self._phi = phi
        self._npix = None
        if type(x)==type(y)==type(z)==type(mu)==type(phi)==np.ndarray:
            assert x.size==y.size==z.size==mu.size==phi.size, 'All input arrays must be of equal size'
            self._npix = x.size
        self._resolution = None
        
    def __getitem__(self, val):
        projection = Projection(
            x=np.array(self._x[val]), 
            y=np.array(self._y[val]),
            z=np.array(self._z[val]), 
            mu=np.array(self._mu[val]), 
            phi=np.array(self._phi[val]),
        )
        return projection
    
    def split(self, n_parts):
        """TODO"""
        x_split = np.array_split(self.x, n_parts) 
        y_split = np.array_split(self.y, n_parts) 
        z_split = np.array_split(self.z, n_parts) 
        mu_split = np.array_split(self.mu, n_parts) 
        phi_split = np.array_split(self.phi, n_parts)        
        projections = [
            Projection(x, y, z, mu, phi) for 
            x, y, z, mu, phi in zip(x_split, y_split, z_split, mu_split, phi_split)
        ]
        return projections
    
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
        return np.rad2deg(np.arccos(self.mu))
 
    @property
    def azimuth(self):
        return np.rad2deg(self.phi)
    
    @property 
    def npix(self):
        return self._npix    
    
    @property
    def resolution(self):
        if self._resolution is None:
            return [self.npix, ]
        else:
            return self._resolution  
    
    
class HomographyProjection(Projection):
    """
    A Homography has a projective tansformation that relates 3D coordinates to pixels.
    """
    def __init__(self):
        super(HomographyProjection, self).__init__()
        
    def project(self, projection_matrix, point_array):
        """
        Project 3D coordinates according to the sensor projection matrix
        
        Parameters
        ----------
        projection matrix: np.array(shape=(3,4), dtype=float)
            The sensor projection matrix K.
        point_array: np.array(shape=(3, num_points), dtype=float)
            An array of num_points 3D points (x,y,z) [km]
        """
        homogenic_point_array = np.pad(point_array,((0,1),(0,0)),'constant', constant_values=1)
        return np.dot(projection_matrix, homogenic_point_array)  
  

class OrthographicProjection(HomographyProjection):
    """
    A parallel ray projection. 
    
    Parameters
    ----------
    bounding_box: shdom.BoundingBox object
        The bounding box is used to compute a projection that will make the entire bounding box visible.
    x_resolution: float
        Pixel resolution [km] in x axis (North)
    y_resolution: float
        Pixel resolution [km] in y axis (East)
    azimuth: float
        Azimuth angle [deg] of the measurements (direciton of the photons)
    zenith: float
        Zenith angle [deg] of the measurements (direciton of the photons)
    altitude: float or 'TOA' (default)
       1. 'TOA': Top of the atmosphere.
       2. float: Altitude of the  measurements.    
    """
    
    def __init__(self, bounding_box, x_resolution, y_resolution, azimuth, zenith, altitude='TOA'):
        super(OrthographicProjection, self).__init__()       
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
    
    
class PerspectiveProjection(HomographyProjection):
    """
    A Perspective trasnormation (pinhole camera).
    
    Parameters
    ----------
    fov: float
        Field of view [deg]
    nx: int
        Number of pixels in camera x axis
    ny: int
        Number of pixels in camera y axis
    x: float
        Location in global x coordinates [km] (North)
    y: float
        Location in global y coordinates [km] (East)
    z: float
        Location in global z coordinates [km] (Up)
    """
    def __init__(self, fov, nx, ny, x, y, z):
        super(PerspectiveProjection, self).__init__()
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
        """
        This is an internal method which is called upon when a rotation matrix is computed to update the global camera coordinates.
        """
        x_c, y_c, z_c = norm(np.matmul(
            self._rotation_matrix, np.matmul(self._inv_k, self._homogeneous_coordinates)))

        self._mu = -z_c.astype(np.float32)
        self._phi = (np.arctan2(y_c, x_c) + np.pi).astype(np.float32)
        self._x = np.full(self.npix, self.position[0], dtype=np.float32)
        self._y = np.full(self.npix, self.position[1], dtype=np.float32)
        self._z = np.full(self.npix, self.position[2], dtype=np.float32)
    
    
    def look_at_transform(self, point, up):
        """
        A look at transform is defined with a point and an up vector.
        
        Parameters
        ----------
        point: np.array(shape=(3,), dtype=float)
            A point in 3D space (x,y,z) coordinates in [km]
        up: np.array(shape=(3,), dtype=float)
            The up vector determines the roll of the camera.
        """
        up = np.array(up)
        direction = np.array(point) - self.position
        zaxis = norm(direction)
        xaxis = norm(np.cross(up, zaxis))
        yaxis = np.cross(zaxis, xaxis)
        self._rotation_matrix = np.stack((xaxis, yaxis, zaxis), axis=1)
        self.update_global_coordinates()
        
        
    def rotate_transform(self, axis, angle):
        """
        Rotate the camera with respect to one of it's (local) axis
        
        Parameters
        ----------
        axis: 'x', 'y' or 'z'
            The rotation axis
        angle: float
            The angle of rotation [deg]
            
        Notes
        -----
        The axis are in the camera coordinates
        """
        
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
    
    
class PrincipalPlaneProjection(Projection):
    """
    Measurments along the principal solar plane.
    
    Parameters
    ----------
    source: shdom.SolarSource
        The source azimuth is used to find the solar principal plane
    x: float
        Location in global x coordinates [km] (North)
    y: float
        Location in global y coordinates [km] (East)
    z: float
        Location in global z coordinates [km] (Up)
    resolution: float 
        Angular resolution of the measurements in [deg]
    """
    def __init__(self, source, x, y, z, resolution=1.0):
        super(PrincipalPlaneProjection, self).__init__()
        self._angles = np.arange(-89.0, 89.0, resolution)
        self._npix = len(self._angles)
        self._x = np.full(self.npix, x, dtype=np.float32)
        self._y = np.full(self.npix, y, dtype=np.float32)
        self._z = np.full(self.npix, z, dtype=np.float32)
        self._mu = (np.cos(np.deg2rad(self._angles))).astype(np.float64)
        self._phi = np.deg2rad(180 * (self._angles < 0.0).astype(np.float) + source.azimuth)
        self._source = source
    
    @property 
    def angles(self):
        return self._angles
        
        
class AlmucantarProjection(Projection):
    """
    Measurments along the solar almucantar.
    
    Parameters
    ----------
    source: shdom.SolarSource
        The source zenith is used to find the solar almucantar plane
    x: float
        Location in global x coordinates [km] (North)
    y: float
        Location in global y coordinates [km] (East)
    z: float
        Location in global z coordinates [km] (Up)
    resolution: float 
        Angular resolution of the measurements in [deg]
    """
    def __init__(self, source, x, y, z, resolution=1.0):
        super(AlmucantarProjection, self).__init__()
        self._phi = np.deg2rad(np.arange(180.0, 360.0, resolution)).astype(np.float64)
        self._npix = len(self._phi)
        self._mu = np.full(self.npix, np.cos(np.deg2rad(source.zenith - 180)), dtype=np.float64)
        self._x = np.full(self.npix, x, dtype=np.float32)
        self._y = np.full(self.npix, y, dtype=np.float32)
        self._z = np.full(self.npix, z, dtype=np.float32)


class HemisphericProjection(Projection):
    """
    Measurments of radiance on a hemisphere.
    
    Parameters
    ----------
    source: shdom.SolarSource
        The source zenith is used to find the solar almucantar plane
    x: float
        Location in global x coordinates [km] (North)
    y: float
        Location in global y coordinates [km] (East)
    z: float
        Location in global z coordinates [km] (Up)
    resolution: float 
        Angular resolution of the measurements in [deg]
    """
    def __init__(self, x, y, z, resolution=5.0):
        super(HemisphericProjection, self).__init__()
        mu = np.cos(np.deg2rad(np.arange(0.0, 80.0+resolution, resolution)))
        phi = np.deg2rad(np.arange(0.0, 360.0+resolution, resolution))
        self._x, self._y, self._z, self._mu, self._phi = np.meshgrid(x, y, z, mu, phi)
        self._x = self._x.ravel().astype(np.float32)
        self._y = self._y.ravel().astype(np.float32)
        self._z = self._z.ravel().astype(np.float32)
        self._mu = self._mu.ravel().astype(np.float64)
        self._phi = self._phi.ravel().astype(np.float64)
        self._npix = self.x.size
        self._resolution = [phi.size, mu.size]
            

class Measurements(object):
    """
    A Measurements object bundles together the imaging geometry and radiance measurents for later optimization.
    It can be initilized with a Camera and images or radiances. 
    Alternatively is can be loaded from file.
    
    Notes
    -----
    When an images exist, radiances are simply a flattened version of the images.
    """
    def __init__(self, camera=None, images=None, radiances=None):
        self._camera = camera
        self._images = images
        self._radiances = radiances
        if images is not None: 
            if type(images) is not list:
                self._images = [images]
            
            #Check if images have the same number of channels
            num_channels_array = np.array(map(lambda img: img.shape[2] if img.ndim>2 else 1, self.images))
            if all([elem == num_channels_array[0] for elem in num_channels_array]):
                num_channels = num_channels_array[0]
            self._radiances = np.concatenate([image.reshape((-1, num_channels), order='F') for image in self.images])
             
    def save(self, path):
        """
        Save Measurements to file.
    
        Parameters
        ----------
        path: str,
            Full path to file. 
        """
        file = open(path,'w')
        file.write(pickle.dumps(self.__dict__, -1))
        file.close()
    
    
    def load(self, path):
        """
        Load Measurements from file.

        Parameters
        ----------
        path: str,
            Full path to file. 
        """        
        file = open(path, 'r')
        data = file.read()
        file.close()
        self.__dict__ = pickle.loads(data)     
    
    def split(self, n_parts):
        """TODO"""
        projections = self.camera.projection.split(n_parts)
        radiances = np.array_split(self.radiances, n_parts) 
        measurements = [shdom.Measurements(
            camera=shdom.Camera(self.camera.sensor, projection), 
            radiances=radiance) for  projection, radiance in zip(projections, radiances)
        ]
        return measurements
    
    def add_noise(self):
        """Add sensor modeled noise to the radiances"""
        raise NotImplemented
    
    @property
    def camera(self):
        return self._camera
    
    @property
    def radiances(self):
        return self._radiances
    
    @property
    def images(self):
        return self._images

    
    
class Camera(object):
    """
    An Camera object ecapsulates both sensor and projection.
    A Sensor needs to define a render method.
    """
    def __init__(self, sensor=Sensor(), projection=Projection()):
        self.set_sensor(sensor)
        self.set_projection(projection)
        
    def set_projection(self, projection):
        self._projection = projection
        
    def set_sensor(self, sensor):
        self._sensor = sensor
        
        # Update function docstring
        if sensor.render.__doc__ is not None:
            self.render.__func__.func_doc += sensor.render.__doc__

    def render(self, rte_solver, n_jobs=1, verbose=0):
        """
        Render an image according to the render function defined by the sensor.
        """
        return self.sensor.render(rte_solver, self.projection, n_jobs, verbose)

    @property
    def projection(self):
        return self._projection
    
    @property
    def sensor(self):
        return self._sensor   

    
    

class MultiViewProjection(Projection):
    """
    A MultiViewProjection object encapsulate several projection geometries for multi-view imaging of a domain.
    
    Parameters
    ----------
    projection_list: list, optional
        A list of Sensor objects
    """
    
    def __init__(self, projection_list=None):
        super(MultiViewProjection, self).__init__()
        self._num_projections = 0
        self._projection_list = []
        self._names = []
        if projection_list:
            for projection in projection_list:
                self.add_projection(projection)

    
    def add_projection(self, projection, name=None):
        """
        Add a projection to the projection list
        
        Parameters
        ----------
        projection: Projection object
            A Projection object to add to the MultiViewProjection
        name: str, optional
            An ID for the projection. 
        """
        attributes = ['x', 'y', 'z', 'mu', 'phi']
        
        if self.num_projections == 0:
            for attr in attributes:
                self.__setattr__('_' + attr, projection.__getattribute__(attr))
            self._npix = [projection.npix]
            self._resolution = [projection.resolution]
            self._names = [name]
        else:
            for attr in attributes:
                self.__setattr__('_' + attr, np.concatenate((self.__getattribute__(attr), 
                                                             projection.__getattribute__(attr))))
            self._npix.append(projection.npix)
            self._names.append(name)  
            self._resolution.append(projection.resolution)
                                    
        self._projection_list.append(projection)
        self._num_projections += 1 

    @property
    def projection_list(self):
        return self._projection_list
    
    @property
    def num_projections(self):
        return self._num_projections    