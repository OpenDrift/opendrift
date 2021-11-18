# Speeding up saving of matplotlib animations to file
#
# Copied from https://stackoverflow.com/questions/30965355/speedup-matplotlib-animation-to-video-file
# Solution by Aule Mahlal - https://stackoverflow.com/users/9291575/aule-mahal
#
# Same method applied also to ImageMagick. Manuel Aghito 18.11.2021
#
# Note: when using iPython this gives an advantage only with qt graphic backend
#       activated by running "%matplotlib qt" does not help with inline graphic backend

from matplotlib.animation import FFMpegWriter
from matplotlib.animation import ImageMagickWriter

class PunkImageMagickWriter(ImageMagickWriter):
    '''Faster ImageMagick-pipe writer bypassing figure.savefig.'''
    def __init__(self, **kwargs):
        '''Initialize the Writer object and sets the default frame_format.'''
        super().__init__(**kwargs)
        self.frame_format = 'rgb' # note rgba doesn't work with ImageMagick

    def grab_frame(self, **savefig_kwargs):
        '''Grab the image information from the figure and save as a movie frame.

        Doesn't use savefig to be faster: savefig_kwargs will be ignored.
        '''
        try:
            # re-adjust the figure size and dpi in case it has been changed by the
            # user.  We must ensure that every frame is the same size or
            # the movie will not save correctly.
            self.fig.set_size_inches(self._w, self._h)
            self.fig.set_dpi(self.dpi)
            # Draw and save the frame as an argb string to the pipe sink
            self.fig.canvas.draw()
            #self._frame_sink().write(self.fig.canvas.tostring_argb()) 
            self._proc.stdin.write(self.fig.canvas.tostring_rgb())
        except (RuntimeError, IOError) as e:
            out, err = self._proc.communicate()
            raise IOError('Error saving animation to file (cause: {0}) '
                      'Stdout: {1} StdError: {2}. It may help to re-run '
                      'with --verbose-debug.'.format(e, out, err))

class PunkFFMpegWriter(FFMpegWriter):
    '''Faster FFMpeg-pipe writer bypassing figure.savefig.'''
    def __init__(self, **kwargs):
        '''Initialize the Writer object and sets the default frame_format.'''
        super().__init__(**kwargs)
        self.frame_format = 'argb'

    def grab_frame(self, **savefig_kwargs):
        '''Grab the image information from the figure and save as a movie frame.

        Doesn't use savefig to be faster: savefig_kwargs will be ignored.
        '''
        try:
            # re-adjust the figure size and dpi in case it has been changed by the
            # user.  We must ensure that every frame is the same size or
            # the movie will not save correctly.
            self.fig.set_size_inches(self._w, self._h)
            self.fig.set_dpi(self.dpi)
            # Draw and save the frame as an argb string to the pipe sink
            self.fig.canvas.draw()
            #self._frame_sink().write(self.fig.canvas.tostring_argb()) 
            self._proc.stdin.write(self.fig.canvas.tostring_argb())
        except (RuntimeError, IOError) as e:
            out, err = self._proc.communicate()
            raise IOError('Error saving animation to file (cause: {0}) '
                      'Stdout: {1} StdError: {2}. It may help to re-run '
                      'with --verbose-debug.'.format(e, out, err))