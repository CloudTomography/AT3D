"""
This module contains summary writer that can be used for monitoring
the progress of long optimizations.
"""
import tensorboardX
import time
import matplotlib.pyplot as plt

class SummaryWriter(tensorboardX.SummaryWriter):
    """
    A wrapper around the tensorboardX.SummaryWriter where any default
    figures/scalars to monitor can be added to avoid defining them
    at the script level.
    """
    def __init__(self, logdir=None, comment='', purge_step=None, max_queue=10, flush_secs=120,
                 filename_suffix='', write_to_disk=True, **kwargs):
        super().__init__(logdir, comment, purge_step, max_queue, flush_secs,
                         filename_suffix, write_to_disk, **kwargs)

    def add_xarray_image(self, tag, image, global_step=None):
        plt.style.use('default')
        image.plot()
        self.add_figure(tag, plt.gcf(), global_step)

    def add_xarray_images(self, tag, image, global_step=None):
        plt.style.use('default')
        image.plot()
        self.add_figure(tag, plt.gcf(), global_step)

    #TODO add some defaults/examples.

class CallbackFn:
    def __init__(self, callback_fn, ckpt_period=-1):
        # callback_fn should return a dictionary.
        self._ckpt_period = ckpt_period
        self._ckpt_time = time.time()
        self._callback_fn = callback_fn
        self.output = {}

    def __call__(self, optimizer=None):
        time_passed = time.time() - self._ckpt_time
        if time_passed > self._ckpt_period:
            self._ckpt_time = time.time()
            out = self._callback_fn(optimizer=optimizer)
            if not self.output:
                for name in out:
                    self.output[name] = []
            for name, value in out.items():
                self.output[name].append(value)
