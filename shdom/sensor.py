"""
Sensor and Sensor related objects used for atmospheric rendering.
The Sensor object defines the projection geometry.
"""

import core
import numpy as np
from shdom import BoundingBox
import itertools
import dill as pickle
from joblib import Parallel, delayed

norm = lambda x: x / np.linalg.norm(x, axis=0)


class Measurements(object):
    """
    A Measurements object bundles together the Sensor geometry and radiance measurents for later optimization.
    TODO
    """
    def __init__(self, sensors=None, **kwargs):
        self._sensors = sensors
        
        if kwargs.has_key('images'):
            self._images = kwargs['images']
            if type(self.images) is list:
                self._radiances = np.concatenate(map(lambda img: img.ravel(order='F'), self.images))
            elif type(self.images) is np.ndarray:
                self._radiances = self.images.ravel(order='F')
        
        if kwargs.has_key('radiances'):
            self._radiances = kwargs['radiances']
        
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
    
    
    def add_noise(self):
        """TODO"""
        raise NotImplemented
    
    
    @property
    def sensors(self):
        return self._sensors
    
    @property
    def radiances(self):
        return self._radiances
    
    @property
    def images(self):
        return self._images
    
    @property
    def wavelength(self):
        return self._wavelength
    
    
class Sensor(object):
    """
    Abstract Sensor class to be inherited by the different types of sensors.
    Each sensor internally defines an arrays of pixel locations (x,y,z) in km and directions (phi, mu).
    TODO
    """
    def __init__(self):
        self._x = None
        self._y = None
        self._z = None
        self._mu = None
        self._phi = None
        self._npix = None
        self.type = 'AbstractSensor'
        self._wavelength = None
     
     
    def render(self, rte_solver):
        """
        The render method integrates a pre-computed in-scatter field (source function) J over the sensor gemoetry.
        The source code for this function is in shdomsub4.f. 
        It is a modified version of the original SHDOM visualize_radiance subroutine in shdomsub2.f.
        
        Parameters
        ----------
        rte_solver: shdom.RteSolver object
            The RteSolver with the precomputed radiative transfer solution (RteSolver.solve method).
        
        Returns
        -------
        radiance: np.array(shape=(sensor.npix,), dtype=np.float)
            The rendered radiances.
        """
        radiance = core.render(
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
            legen=rte_solver._legen.reshape(rte_solver._nleg+1, -1),            
            dirflux=rte_solver._dirflux,
            fluxes=rte_solver._fluxes,
            source=rte_solver._source,          
            srctype=rte_solver._srctype,
            sfctype=rte_solver._sfctype,
            units=rte_solver._units
        )
        return radiance


    def par_render(self, rte_solver, n_jobs=30, verbose=0):
        """
        The par_render method integrates a pre-computed in-scatter field (source function) J over the sensor gemoetry.
        The source code for this function is in shdomsub4.f. 
        It is a modified version of the original SHDOM visualize_radiance subroutine in shdomsub2.f.
        This method is a parallel version of the render method.
        
        Parameters
        ----------
        rte_solver: shdom.RteSolver object
            The RteSolver with the precomputed radiative transfer solution (RteSolver.solve method).
        n_jobs: int, default=30
            The number of jobs to divide the rendering into.
        verbose: int, default=0
            How much verbosity in the parallel rendering proccess.
            
        Returns
        -------
        radiance: np.array(shape=(sensor.npix,), dtype=np.float)
            The rendered radiances.
        
        Notes
        -----
        For a small domain, or small amout of pixels, par_render is slower than render.
        """
        maxnlm=16384
        maxleg=2000
        maxphase=500000
        maxscatang=361
        nscatangle = max(36, min(maxscatang, 2*rte_solver._nleg))
        assert ((rte_solver._ml+1)**2-(2-rte_solver._ncs)*(rte_solver._ml*(rte_solver._ml+1))/2 < maxnlm), '[par_visualize1] assert: maxnlm exceeded'        
        assert (rte_solver._nleg < maxleg), '[par_render] assert: maxleg exceeded' 
        assert (rte_solver._pa.numphase < maxphase), '[par_render] assert: maxphase exceeded' 
        
        if rte_solver._srctype is not 'T':
            rte_solver._ylmsun = core.ylmall(
                mu=rte_solver._solarmu,
                phi=rte_solver._solaraz,
                ml=rte_solver._ml,
                mm=rte_solver._mm,
                ncs=rte_solver._ncs,
                p=rte_solver._ylmsun
            )
            
            if rte_solver._deltam and rte_solver._pa.numphase > 0:
                rte_solver._phasetab = core.precompute_phase(
                    maxphase=maxphase,
                    nscatangle=nscatangle,
                    numphase=rte_solver._pa.numphase,
                    ml=rte_solver._ml,
                    nleg=rte_solver._nleg,
                    legen=rte_solver._legen.reshape(rte_solver._nleg+1, -1)
                )  
            
        # Make the isotropic radiances for the top boundary
        rte_solver._bcrad = core.compute_top_radiances(
            srctype=rte_solver._srctype,
            skyrad=rte_solver._skyrad,
            waveno=rte_solver._waveno,
            wavelen=rte_solver._wavelen,
            units=rte_solver._units,
            ntoppts=rte_solver._ntoppts,
            bcrad=rte_solver._bcrad
        )
        
        # Make the bottom boundary radiances for the Lambertian surfaces.  
        # Compute the upwelling bottom radiances using the downwelling fluxes.        
        if rte_solver._sfctype == 'FL':
            rte_solver._bcrad = core.fixed_lambertian_boundary(
                nbotpts=rte_solver._nbotpts,
                bcptr=rte_solver._bcptr,
                dirflux=rte_solver._dirflux,
                fluxes=rte_solver._fluxes,
                srctype=rte_solver._srctype,
                gndtemp=rte_solver._gndtemp,
                gndalbedo=rte_solver._gndalbedo,
                waveno=rte_solver._waveno,
                wavelen=rte_solver._wavelen,
                units=rte_solver._units,
                bcrad=rte_solver._bcrad
            )
            
        elif rte_solver._sfctype == 'VL':
            raise NotImplementedError('Variable surface not implemented.')       
        
        x_split = np.array_split(self.x, n_jobs) 
        y_split = np.array_split(self.y, n_jobs) 
        z_split = np.array_split(self.z, n_jobs) 
        mu_split = np.array_split(self.mu, n_jobs) 
        phi_split = np.array_split(self.phi, n_jobs)        
        npix_split = map(lambda x: len(x), x_split)
        
        radiance = Parallel(n_jobs=n_jobs, backend="threading", verbose=verbose)(
            delayed(core.par_render, check_pickle=False)(
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
                legen=rte_solver._legen.reshape(rte_solver._nleg+1, -1),            
                dirflux=rte_solver._dirflux,
                fluxes=rte_solver._fluxes,
                source=rte_solver._source, 
                camx=x,
                camy=y,
                camz=z,
                cammu=mu,
                camphi=phi, 
                npix=npix,
                srctype=rte_solver._srctype,
                sfctype=rte_solver._sfctype,
                units=rte_solver._units,
                phasetab=rte_solver._phasetab,
                ylmsun=rte_solver._ylmsun,
                nscatangle=nscatangle
            ) for x, y, z, mu, phi, npix in zip(x_split, y_split, z_split, mu_split, phi_split, npix_split))   
        
        return np.concatenate(radiance)
    
    
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
    def npix(self):
        return self._npix    
    
    @property
    def wavelength(self):
        return self._wavelength
    
    
class ProjectiveMonochromeSensor(Sensor):
    """
    A ProjectiveMonochromeSensor has a projective tansformation that relates 3D coordinate to pixels.
    When rendering, this sensor will produce 2D images.
    TODO
    """
    def __init__(self, wavelength):
        super(ProjectiveMonochromeSensor, self).__init__()
        self._type = 'ProjectiveMonochromeSensor'
        self._resolution = None
        self._wavelength = wavelength
    
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
        
    def render(self, rte_solver):
        """
        The render method integrates a pre-computed in-scatter field (source function) J over the sensor gemoetry.
        The source code for this function is in shdomsub4.f. 
        It is a modified version of the original SHDOM visualize_radiance subroutine in shdomsub2.f.
        
        Parameters
        ----------
        rte_solver: shdom.RteSolver object
            The RteSolver with the precomputed radiative transfer solution (RteSolver.solve method).
        
        Returns
        -------
        radiance: np.array(shape=(sensor.resolution), dtype=np.float)
            The rendered radiances.
        """
        visout = super(ProjectiveMonochromeSensor, self).render(rte_solver)
        return visout.reshape(self.resolution, order='F')
    
    def par_render(self, rte_solver, n_jobs=30, verbose=0):
        """
        The par_render method integrates a pre-computed in-scatter field (source function) J over the sensor gemoetry.
        The source code for this function is in shdomsub4.f. 
        It is a modified version of the original SHDOM visualize_radiance subroutine in shdomsub2.f.
        This method is a parallel version of the render method.
        
        Parameters
        ----------
        rte_solver: shdom.RteSolver object
            The RteSolver with the precomputed radiative transfer solution (RteSolver.solve method).
        n_jobs: int, default=30
            The number of jobs to divide the rendering into.
        verbose: int, default=0
            How much verbosity in the parallel rendering proccess.
            
        Returns
        -------
        radiance: np.array(shape=(sensor.resolution), dtype=np.float)
            The rendered radiances.
        
        Notes
        -----
        For a small domain, or small amout of pixels, par_render is slower than render.
        """        
        visout = super(ProjectiveMonochromeSensor, self).render(rte_solver, n_jobs, verbose)
        return np.array(visout).reshape(self.resolution, order='F')
    
    
    @property
    def resolution(self):
        return self._resolution
    
    
class SensorArray(Sensor):
    """
    A SensorArray object encapsulate several sensors e.g. for multi-view imaging of a domain.
    
    Parameters
    ----------
    sensor_list: list, optional
        A list of Sensor objects
    """
    
    def __init__(self, sensor_list=None):
        super(SensorArray, self).__init__()
        self.type = 'SensorArray'
        self._num_sensors = 0
        self._sensor_list = []
        if sensor_list:
            for sensor in sensor_list:
                self.add_sensor(sensor)
    
    
    def add_sensor(self, sensor, id=None):
        """
        Add a sensor to the SensorArray
        
        Parameters
        ----------
        sensor: Sensor object
            A Sensor object to add to the SensorArray
        id: str, optional
            An ID for the sensor. 
        """
        attributes = ['x', 'y', 'z', 'mu', 'phi']
        
        if self.num_sensors is 0:
            for attr in attributes:
                self.__setattr__('_' + attr, sensor.__getattribute__(attr))
            self._npix = [sensor.npix]
            self._resolution = [sensor.resolution]
            self._type = [sensor.type]
            self._ids = ['sensor0' if id is None else id]
        else:
            for attr in attributes:
                self.__setattr__('_' + attr, np.concatenate((self.__getattribute__(attr), 
                                                             sensor.__getattribute__(attr))))
            self._npix.append(sensor.npix)
            self._resolution.append(sensor.resolution)
            self._type.append(sensor.type)
            self._ids.append('sensor'+str(self.num_sensors) if id is None else id)  
        
        self._sensor_list.append(sensor)
        self._num_sensors += 1
    
    
    def render(self, rte_solver):
        """
        Serial rendering of each sensor.
        
        Parameters
        ----------
        rte_solver: shdom.RteSolver object
            The RteSolver with the precomputed radiative transfer solution (RteSolver.solve method).
        
        Returns
        -------
        radiance: list
            A list of the rendered radiances per sensor.
        """
        split_measurements = [sensor.render(rte_solver) for sensor in self.sensor_list]
        radiance = []
        for sensor, meas in zip(self.sensor_list, split_measurements):
            meas = np.array(meas)
            if hasattr(sensor, 'resolution'):
                meas = meas.reshape(sensor.resolution, order='F')
            radiance.append(meas)            
        return radiance     


    def par_render(self, rte_solver, n_jobs=30, verbose=0):
        """
        Parallel rendering of all sensors.
        
        Parameters
        ----------
        rte_solver: shdom.RteSolver object
            The RteSolver with the precomputed radiative transfer solution (RteSolver.solve method).
        
        Returns
        -------
        radiance: list
            A list of the rendered radiances per sensor.
        """
        radiance = super(SensorArray, self).par_render(rte_solver, n_jobs, verbose)
        
        split_measurements = np.split(radiance, np.cumsum(self.npix[:-1]))
        radiance = []
        for sensor, meas in zip(self.sensor_list, split_measurements):
            meas = np.array(meas)
            if hasattr(sensor, 'resolution'):
                meas = meas.reshape(sensor.resolution, order='F')
            radiance.append(meas)
        return radiance    


    @property
    def sensor_list(self):
        return self._sensor_list
    
    @property
    def num_sensors(self):
        return self._num_sensors    
    
    
class OrthographicMonochromeSensor(ProjectiveMonochromeSensor):
    """
    A parallel ray projection Sensor measures radiance along rays. 
    
    Parameters
    ----------
    wavelength: float
        Monochromatic wavelength [Micron]
    bounding_box: shdom.BoundingBox object
        The bounding box is used to compute a projection that will make the entire bounding box visible.
    x_resolution: float
        Pixel resolution [km] in x axis (North)
    y_resolution: float
        Pixel resolution [km] in y axis (East)
    azimuth: float
        Azimuth angle [deg] of the radiance measurements (direciton of the photons)
    zenith: float
        Zenith angle [deg] of the radiance measurements (direciton of the photons)
    altitude: float or 'TOA' (default)
       1. 'TOA': measurements of exiting radiace at the top of the atmosphere.
       2. float: Altitude of the radiance measurements.    
    """
    
    def __init__(self, wavelength, bounding_box, x_resolution, y_resolution, azimuth, zenith, altitude='TOA'):
        super(OrthographicMonochromeSensor, self).__init__(wavelength)
        self.type = 'OrthographicMonochromeSensor'
        
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
    
    

class PerspectiveMonochromeSensor(ProjectiveMonochromeSensor):
    """
    A Monochromatic Sensor with a Perspective trasnormation 
    (a pinhole camera).
    
    wavelength: float
        Monochromatic wavelength [Micron]
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
    
    def __init__(self, wavelength, fov, nx, ny, x, y, z):
        super(PerspectiveMonochromeSensor, self).__init__(wavelength)
        self.type = 'PerspectiveMonochromeSensor'

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
