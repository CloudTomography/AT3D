"""
Visualization tools
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interact
from mpl_toolkits.axes_grid1 import make_axes_locatable

@xr.register_dataarray_accessor("visualization")
class _VisualizationAccessor(object):
    """
    Register a custom accessor VisualizationAccessor on xarray.DataArray object.
    This adds methods for visualization of scalar valued DataArrays.
    """
    def __init__(self, xarray_obj):
        self._obj = xarray_obj

    def slider(self, dim, ax=None, cmap=None):
        """
        Interactive slider visualization of a 3D xr.DataArray along specified dimension.

        Parameters
        ----------
        dim: str,
            The dimension along which to visualize frames
        ax: matplotlib axis,
            A matplotlib axis object for the visualization.
        cmap : str or matplotlib.colors.Colormap, optional
            The Colormap instance or registered colormap name used to map scalar data to colors.
            Defaults to :rc:`image.cmap`.
        """
        data = self._obj.squeeze()
        if data.ndim != 3:
            raise AttributeError('Move dimensions ({}) different than 3'.format(data.ndim))

        num_frames = data[dim].size
        image_dims = list(data.dims)
        image_dims.remove(dim)

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()

        extent = [data[image_dims[0]].min(), data[image_dims[0]].max(),
                  data[image_dims[1]].min(), data[image_dims[1]].max()]

        im = ax.imshow(data.isel({dim: 0}), extent=extent, cmap=cmap)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cbar = fig.colorbar(im, cax=cax)

        def imshow_frame(frame):
            img = data.isel({dim: frame})
            ax.imshow(img, origin='lower', extent=extent, cmap=cmap)
            cbar.mappable.set_clim([img.min(), img.max()])
            ax.set_title('{} = {:2.2f}'.format(dim, float(img[dim])))

        interact(imshow_frame, frame=(0, num_frames-1));
