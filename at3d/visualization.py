"""
Visualization tools
"""
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interact
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import animation

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
            ax.imshow(img.T, origin='lower', extent=extent, cmap=cmap)
            cbar.mappable.set_clim([img.min(), img.max()])
            ax.set_title('{} = {:2.2f}'.format(dim, float(img[dim])))

        interact(imshow_frame, frame=(0, num_frames-1));

    def animate(self, dim, ax=None, vmin=None, vmax=None, cmap='RdBu_r', add_ticks=True, add_colorbar=True,
                fps=10, output=None):
        """
        Animate a 3D xr.DataArray along a chosen dimension.

        Parameters
        ----------
        dim: str,
            The dimension along which to animate frames
        ax: matplotlib axis,
            A matplotlib axis object for the visualization.
        vmin, vmax : float, optional
            vmin and vmax define the data range that the colormap covers.
            By default, the colormap covers the complete value range of the supplied data.
        cmap : str or matplotlib.colors.Colormap, default='RdBu_r'
            The Colormap instance or registered colormap name used to map scalar data to colors.
            Defaults to :rc:`image.cmap`.
        add_ticks: bool, default=True
            If true then ticks will be visualized.
        add_colorbar: bool, default=True
            If true then a colorbar will be visualized
        fps: float, default=10,
            Frames per seconds.
        output: string,
            Path to save the animated gif. Should end with .gif.

        Returns
        -------
        anim: matplotlib.animation.FuncAnimation
            Animation object.
        """
        movie = self._obj.squeeze()
        if movie.ndim != 3:
            raise AttributeError('Movie dimensions ({}) different than 3'.format(movie.ndim))

        num_frames = movie[dim].size
        image_dims = list(movie.dims)
        image_dims.remove(dim)
        nx, ny = [movie.sizes[dim] for dim in image_dims]

        # Image animation function (called sequentially)
        def animate_frame(i):
            im.set_array(movie.isel({dim: i}).T)
            ax.set_title('{} = {:2.2f}'.format(dim, float(movie[dim].isel({dim: i}))))
            return [im]

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()

        extent = [movie[image_dims[0]].min(), movie[image_dims[0]].max(),
                  movie[image_dims[1]].min(), movie[image_dims[1]].max()]

        # Initialization function: plot the background of each frame
        im = ax.imshow(np.zeros((nx, ny)), extent=extent, origin='lower', cmap=cmap)
        if add_colorbar:
            fig.colorbar(im)
        if add_ticks == False:
            ax.set_xticks([])
            ax.set_yticks([])

        if (vmin is None) & (vmax is None):
            vmin = movie.min()
            vmax = movie.max()
            if vmin < 0:
                vmin = -1*max(np.abs(vmax), np.abs(vmin))
                vmax = max(np.abs(vmax), np.abs(vmin))
        elif vmin is None:
            vmin = movie.min()
        elif vmax is None:
            vmax = movie.max()

        im.set_clim(vmin, vmax)
        anim = animation.FuncAnimation(fig, animate_frame, frames=num_frames, interval=1e3 / fps)

        if output is not None:
            anim.save(output, writer='imagemagick', fps=fps)
        return anim
