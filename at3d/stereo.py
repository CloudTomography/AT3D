"""
This module contains an interface to the More Global Matching (MGM) software
that is publicly available. It does not implement stereo matching itself.
You need to install the MGM matcher first from:

https://github.com/gfacciol/mgm

If the MGM package doesn't work you are on your own.

Cite the software if you use it:

@inproceedings{BMVC2015_90,
	title={MGM: A Significantly More Global Matching for Stereovision},
	author={Gabriele Facciolo and Carlo de Franchis and Enric Meinhardt},
	year={2015},
	month={September},
	pages={90.1-90.12},
	articleno={90},
	numpages={12},
	booktitle={Proceedings of the British Machine Vision Conference (BMVC)},
	publisher={BMVA Press},
	editor={Xianghua Xie, Mark W. Jones, and Gary K. L. Tam},
	doi={10.5244/C.29.90},
	isbn={1-901725-53-7},
	url={https://dx.doi.org/10.5244/C.29.90}
}
"""

import os
import datetime
import subprocess
import sys
import shutil
import warnings

import numpy as np
import xarray as xr
from PIL import Image


def digitize_image(img, dtype=np.uint8):
    """Round an image to the specified digitization accuracy.

    Parameters
    ----------
    img : np.ndarray
        The data to be digitized into an integer precision.
    dtype : type
        The numerical type of the integer precision to digitize the image to.

    Returns
    -------
    rounded_im : np.ndarray
        The digitized output.
    """
    p = 0
    a, b = np.nanpercentile(img, (p, 100 - p))
    rounded_im = np.round(np.iinfo(dtype).max * (np.clip(img, a, b) - a) / (b - a)).astype(dtype)
    return rounded_im

# This function is adapted from https://github.com/centreborelli/s2p
def run(cmd, env=os.environ, timeout=None, shell=False, verbose=False):
    """
    Runs a shell command, and print it before running.
    Arguments:
        cmd: list of a command and its arguments, or as a fallback,
            a string to be passed to a shell that will be split into a list.
        env (optional, default value is os.environ): dictionary containing the
            environment variables
        timeout (optional, int): time in seconds after which the function will
            raise an error if the command hasn't returned
        TODO: remove the temporary `shell` argument once all commands use shell=False
        shell (bool): run the command in a subshell. Defaults to False.
    Both stdout and stderr of the shell in which the command is run are those
    of the parent process.
    """
    if verbose:
        print("\nRUN: %s" % cmd)
    t = datetime.datetime.now()
    if not isinstance(cmd, list) and not shell:
        cmd = cmd.split()
    subprocess.run(cmd, shell=shell, stdout=sys.stdout, stderr=sys.stderr,
                   env=env, timeout=timeout, check=True)
    if verbose:
        print(datetime.datetime.now() - t)

class StereoMatcher:
    """
    Initializes a "More Global Matching" (MGM) matcher that can be used to
    calculate correspondences between two different images.

    Parameters
    ----------
    mgm_directory : str
        The path to the directory in which the "mgm" executable is contained.
    temp_directory : str
        The name of a temporary directory which will be created inside `mgm_directory`
        and used to store output from the MGM matcher.
    matching_cost : str
        The matching cost function e.g. normalized cross correlation or
        census transform.
    P1 : int
        A regularization term that penalizes having a difference in disparity between
        adjacent pixels that is larger than 1.
    P2 : int
        A regularization term that penalizes having a difference in disparity between
        adjacent pixels that is larger than 2.
    n_directions : int
        The number of 'semi-global' paths to use to calculate the regularization term.
        More paths results in a better regularization term.
    subpixel : string
        The type of sub-pixel refinement of the disparity.
    OMP_NUM_THREADS : int
        The number of threads to use. Set as an environmental variables
    MEDIAN : int
        The size of a median filter to post-process the disparities.
        Set as an environmental variable.
    CENSUS_NCC_WIN : int
        The window size for the matching cost in pixels.
        Set as an environmental variable.
    TSGM : int
        The level of regularity. See MGM paper for more details.
        Set as an environmental variable.
	prefilter : str
		'none', 'census', 'sobelx', 'gblur'
    """
    def __init__(self, mgm_directory='./', temp_directory='./mgm_temp',
                matching_cost='census', P1=1, P2=10,
                n_directions=8, subpixel='cubic',
                OMP_NUM_THREADS=4, MEDIAN=1,
                CENSUS_NCC_WIN=5, TSGM=3, prefilter='none',
				aP1=1, aP2=1, aThresh=5):

        # do some checks that mgm is there.
        joined_path = os.path.join(mgm_directory, 'mgm')
        if (not 'mgm' in os.listdir(mgm_directory)) or  \
            (shutil.which(joined_path) != joined_path):
            raise OSError(
            "The mgm executable was not found in the directory '{}'".format(mgm_directory)
            )

        self._mgm_directory = mgm_directory

        # where temporary images will be stored.
        if not os.path.exists(os.path.join(os.getcwd(), temp_directory)):
            os.makedirs(os.path.join(os.getcwd(), temp_directory))
        self._temp_directory = os.path.join(os.getcwd(), temp_directory)

        config = {
            'OMP_NUM_THREADS': str(OMP_NUM_THREADS),
            'MEDIAN': str(MEDIAN),
            'CENSUS_NCC_WIN': str(CENSUS_NCC_WIN),
            'TSGM': str(TSGM)
            }
        self._config = config
        self.config = config
        self.matching_cost = matching_cost
        self.P1 = P1
        self.P2 = P2
        self.n_directions = n_directions
        self.subpixel = subpixel
        self.prefilter = prefilter
        self.aP1 = aP1
        self.aP2 = aP2
        self.aThresh = aThresh

    @property
    def config(self):
        """View the current configuration of the matcher."""
        return self._config

    @config.setter
    def config(self, config):
        """Update the configuration of the matcher."""
        for var, val in config.items():
            os.environ[var] = val
            self._config[var] = val

    def match_stereorectified_images(self, ref_image, other_image, min_disparity=-20,
									 max_disparity=20, verbose=False):
        """
        Run the matcher on two images.

        Runs the MGM matcher as a subprocess. Matching occurs
        along the second dimension of the images. Images are assumed to already be epipolar
        rectified. The search window is one dimensional. For more details read
        the documentation for the mgm matcher.

        Parameters
        ----------
        ref_image : xr.Dataset, xr.DataArray or np.ndarray
            This is the image that disparities are calculated relative to.
            Either a dataset with a `I` variable that is used 2D and used as the image
            or a numpy array.
        other_image :
            This is the comparison image.
            Either a dataset with a `I` variable that is used 2D and used as the image
            or a numpy array.
	    min_disparity : int
	        The minimum disparity to search for valid correspondences between images
	        in units of pixels in the image.
	    max_disparity : int
	        The maximum disparity to search for valid correspondences between images
	        in units of pixels in the image.

        Returns
        -------
        disparity : np.ndarray
            The 2D image (same as `ref_image`) which has the pixel disparities
            to the corresponding best matching points in `other_image`.
        cost : np.ndarray
            The cost of these best-fitting disparities.
        backflow : np.ndarray
            The values of `other_image` backflowed to the `ref_image`.
        """
        current_cwd = os.getcwd()
        os.chdir(self._mgm_directory)

        ref_other = []
        for name, image in zip(['ref_image', 'other_image'], [ref_image, other_image]):
            if isinstance(image, xr.Dataset):
                img = image.I.data
            elif isinstance(image, xr.DataArray):
                img = image.data
            elif isinstance(image, np.ndarray):
                img = image
            else:
                raise TypeError(
                    """`{}` should be a numpy array or xr.Dataset/xr.DataArray""".format(name)
                )

            ref_other.append(digitize_image(img, dtype=np.uint16))

        full_names = []
        for name, im in zip(['ref_image', 'other_image'], ref_other):

            im = Image.fromarray(im)
            full_name = os.path.join(self._temp_directory, name + '.tif')
            im.save(full_name)
            full_names.append(full_name)

        disparity_name = os.path.join(self._temp_directory, 'disparity.tif')
        cost_name = os.path.join(self._temp_directory, 'cost.tif')
        backflow_name = os.path.join(self._temp_directory, 'backflow.tif')

        cmd = """./mgm -r {0} -R {1} -s {2} -t {3} -O {4} -P1 {5} -P2 {6} -p {7} -aP1 {8} -aP2 {9} -aThresh {10} """.format(
            min_disparity, max_disparity, self.subpixel, self.matching_cost,
            self.n_directions, self.P1, self.P2, self.prefilter, self.aP1, self.aP2, self.aThresh \
        ) + full_names[0] + ' ' + full_names[1] + ' ' + disparity_name + \
        ' ' + cost_name +' ' + backflow_name

        run(cmd, verbose=verbose)

        disparity = np.array(Image.open(disparity_name))
        cost = np.array(Image.open(cost_name))
        backflow = np.array(Image.open(backflow_name))

		# undo digitization on backflow
        a, b = np.nanpercentile(other_image, (0, 100))
        backflow = (backflow*(b-a)/np.iinfo(np.uint16).max) + a

        os.remove(disparity_name)
        os.remove(cost_name)
        os.remove(backflow_name)
        os.remove(full_names[0])
        os.remove(full_names[1])

        os.chdir(current_cwd)

        return disparity, cost, backflow
