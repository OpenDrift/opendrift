import os
import logging; logging.captureWarnings(True)
import warnings

import numpy as np

import matplotlib
matplotlib.rcParams['legend.numpoints'] = 1
matplotlib.rcParams['legend.scatterpoints'] = 1
import matplotlib.pyplot as plt
from matplotlib import animation
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

class Plotter:
    def set_up_map(self, corners=None, buffer=.1, delta_lat=None,
                   lscale=None, fast=False, hide_landmask=False, **kwargs):
        """
        Generate Figure instance on which trajectories are plotted.

        :param hide_landmask: do not plot landmask (default False)
        :type hide_landmask: bool

        provide corners=[lonmin, lonmax, latmin, latmax] for specific map selection
        """

        # Initialise map
        if corners is not None:
            lons, lats = self.get_lonlats()  # TODO: to be removed
            lonmin = corners[0]
            lonmax = corners[1]
            latmin = corners[2]
            latmax = corners[3]
        elif hasattr(self, 'ds'):
            lons = self.ds.lon
            lats = self.ds.lat
            if hasattr(self, 'lonmin'):
                lonmin = self.lonmin
                lonmax = self.lonmax
                latmin = self.latmin
                latmax = self.latmax
            else:
                self.logger.debug('Finding min longitude...')
                lonmin = np.nanmin(self.ds.lon)
                self.logger.debug('Finding max longitude...')
                lonmax = np.nanmax(self.ds.lon)
                self.logger.debug('Finding min latitude...')
                latmin = np.nanmin(self.ds.lat)
                self.logger.debug('Finding max latitude...')
                latmax = np.nanmax(self.ds.lat)
                self.lonmin = lonmin
                self.lonmax = lonmax
                self.latmin = latmin
                self.latmax = latmax
                if os.path.exists(self.analysis_file):
                    self.af = Dataset(self.analysis_file, 'a')
                else:
                    self.af = Dataset(self.analysis_file, 'w')
                self.af.lonmin = lonmin
                self.af.lonmax = lonmax
                self.af.latmin = latmin
                self.af.latmax = latmax
                self.af.close()
        else:
            lons, lats = self.get_lonlats()
            lonmin = np.nanmin(lons) - buffer*2
            lonmax = np.nanmax(lons) + buffer*2
            latmin = np.nanmin(lats) - buffer
            latmax = np.nanmax(lats) + buffer

        if fast is True:
            self.logger.warning("plotting fast. This will make your plots less accurate.")

            import matplotlib.style as mplstyle
            mplstyle.use(['fast'])

            # use a spherical earth
            axis = 57.29577951308232  # something to do with pi
            globe = ccrs.Globe(ellipse = None, semimajor_axis = axis, semiminor_axis = axis)
            crs = ccrs.Mercator(globe = globe)

            if lscale is None:
                lscale = 'c'
        else:
            crs = ccrs.Mercator()
            if lscale is None:
                lscale = 'auto'

        meanlat = (latmin + latmax)/2
        aspect_ratio = np.float(latmax-latmin) / (np.float(lonmax-lonmin))
        aspect_ratio = aspect_ratio / np.cos(np.radians(meanlat))
        if aspect_ratio > 1:
            fig = plt.figure(figsize=(11./aspect_ratio, 11.))
        else:
            fig = plt.figure(figsize=(11., 11.*aspect_ratio))

        ax = fig.add_subplot(111, projection=crs)  # need '111' for Python 2
        ax.set_extent([lonmin, lonmax, latmin, latmax], crs=ccrs.PlateCarree())

        def show_landmask(landmask):
            maxn = 512.
            dx = (lonmax - lonmin) / maxn
            dy = (latmax - latmin) / maxn
            dx = max(landmask.dx, dx)
            dy = max(landmask.dy, dy)

            x = np.array([lonmin, lonmax])
            y = np.array([latmin, latmax])
            ndx = np.int32(dx / landmask.dx)  # normalized delta
            ndy = np.int32(dy / landmask.dy)

            xm, ym = landmask.invtransform * (x, y)
            xm = xm.astype(np.int32)
            ym = ym.astype(np.int32)

            img = landmask.mask[ym[0]:ym[1]:ndy, xm[0]:xm[1]:ndx]

            from matplotlib import colors
            cmap = colors.ListedColormap(['white', 'gray'])
            ax.imshow(img, origin = 'lower', extent=[lonmin, lonmax, latmin, latmax],
                      transform=ccrs.PlateCarree(), cmap=cmap)

        import shapely
        from shapely.geometry import box

        if not hide_landmask:
            if 'land_binary_mask' in self.priority_list and (
                self.priority_list['land_binary_mask'][0] == 'shape' or (
                self.priority_list['land_binary_mask'][0] == 'global_landmask' \
                and not self.readers['global_landmask'].skippoly \
                and (self.readers['global_landmask'].mask.extent is None \
                    or self.readers['global_landmask'].mask.extent.contains(box(lonmin, latmin, lonmax, latmax))))):

                if self.priority_list['land_binary_mask'][0] == 'global_landmask':
                    self.logger.debug("Using existing GSHHS shapes..")
                    landmask = self.readers['global_landmask'].mask
                    if fast:
                        show_landmask(landmask)

                    else:
                        extent = box(lonmin, latmin, lonmax, latmax)
                        extent = shapely.prepared.prep(extent)
                        polys = [p for p in landmask.polys.geoms if extent.intersects(p)]

                        ax.add_geometries(polys,
                                ccrs.PlateCarree(),
                                facecolor=cfeature.COLORS['land'],
                                edgecolor='black')
                else:  # Using custom shape reader
                    ax.add_geometries(self.readers['shape'].polys, ccrs.PlateCarree(), facecolor=cfeature.COLORS['land'], edgecolor='black')

            else:
                if fast:
                    from opendrift_landmask_data import Landmask
                    show_landmask(Landmask(skippoly=True))

                else:
                    self.logger.debug ("Adding GSHHS shapes..")
                    f = cfeature.GSHHSFeature(scale=lscale, levels=[1],
                            facecolor=cfeature.COLORS['land'])
                    ax.add_geometries(
                            f.intersecting_geometries([lonmin, lonmax, latmin, latmax]),
                            ccrs.PlateCarree(),
                            facecolor=cfeature.COLORS['land'],
                            edgecolor='black')


        gl = ax.gridlines(ccrs.PlateCarree(), draw_labels=True)
        if cartopy.__version__ < '0.18.0':
            gl.xlabels_top = False  # Cartopy < 0.18
        else:
            gl.top_labels = None  # Cartopy >= 0.18

        fig.canvas.draw()
        fig.set_tight_layout(True)

        if not hasattr(self, 'ds'):
            try:
                firstlast = np.ma.notmasked_edges(lons, axis=1)
                index_of_first = firstlast[0][1]
                index_of_last = firstlast[1][1]
            except:
                index_of_last = 0
        else:
            index_of_first = None
            index_of_last = None

        try:  # Activate figure zooming
            mng = plt.get_current_fig_manager()
            mng.toolbar.zoom()
        except:
            pass

        try:  # Maximise figure window size
            mng.resize(*mng.window.maxsize())
        except:
            pass

        return fig, ax, crs, lons.T, lats.T, index_of_first, index_of_last
    def _get_comparison_xy_for_plots(self, compare):
        if not type(compare) is list:
            compare = [compare]
        compare_list = [{}]*len(compare)
        for cn, comp in enumerate(compare):
            compare_list[cn] = {}
            cd = compare_list[cn]  # pointer to dict with data
            if type(comp) is str:
                # Other is given as filename
                other = self.__class__(loglevel=0)
                other.io_import_file(comp)
            else:
                # Other is given as an OpenDrift object
                other = comp

            # Find map coordinates of comparison simulations
            cd['x_other'], cd['y_other'] = \
                (other.history['lon'].copy(), other.history['lat'].copy())
            cd['x_other_deactive'], cd['y_other_deactive'] = \
                (other.elements_deactivated.lon.copy(),
                    other.elements_deactivated.lat.copy())
            cd['firstlast'] = np.ma.notmasked_edges(
                cd['x_other'], axis=1)
            cd['index_of_last_other'] = cd['firstlast'][1][1]
            cd['index_of_last_deactivated_other'] = \
                cd['index_of_last_other'][other.elements_deactivated.ID-1]

        return compare_list

    def plot(self, background=None, buffer=.2, corners=None, linecolor=None, filename=None,
             show=True, vmin=None, vmax=None, compare=None, cmap='jet',
             lvmin=None, lvmax=None, skip=2, scale=10, show_scalar=True,
             contourlines=False, trajectory_dict=None, colorbar=True,
             linewidth=1, lcs=None, show_particles=True, show_initial=True,
             density_pixelsize_m=1000,
             surface_color=None, submerged_color=None, markersize=20,
             title='auto', legend=True, legend_loc='best', lscale=None,
             fast=False, hide_landmask=False, **kwargs):
        """Basic built-in plotting function intended for developing/debugging.

        Plots trajectories of all particles.
        Positions marked with colored stars:
        - green: all start positions
        - red: deactivated particles
        - blue: particles still active at end of simulation

        Requires availability of Cartopy.

        Arguments:
            background: string, name of variable (standard_name) which will
                be plotted as background of trajectories, provided that it
                can be read with one of the available readers.

            buffer: float; spatial buffer of plot in degrees of
                longitude/latitude around particle collection.

            background: name of variable to be plotted as background field.
            Use two element list for vector fields, e.g. ['x_wind', 'y_wind']

            vmin, vmax: minimum and maximum values for colors of background.

            linecolor: name of variable to be used for coloring trajectories, or matplotlib color string.

            lvmin, lvmax: minimum and maximum values for colors of trajectories.

            lscale (string): resolution of land feature ('c', 'l', 'i', 'h', 'f', 'auto'). default is 'auto'.

            fast (bool): use some optimizations to speed up plotting at the cost of accuracy

            :param hide_landmask: do not plot landmask (default False). See :ref:`model_landmask_only_model` for example usage.
            :type hide_landmask: bool
        """


        mappable = None

        if self.num_elements_total() == 0:
            raise ValueError('Please run simulation before plotting')

        start_time = datetime.now()

        # x, y are longitude, latitude -> i.e. in a PlateCarree CRS
        gcrs = ccrs.PlateCarree()
        fig, ax, crs, x, y, index_of_first, index_of_last = \
            self.set_up_map(buffer=buffer, corners=corners, lscale=lscale, fast=fast, hide_landmask=hide_landmask, **kwargs)

        markercolor = self.plot_comparison_colors[0]

        # The more elements, the more transparent we make the lines
        min_alpha = 0.1
        max_elements = 5000.0
        alpha = min_alpha**(2*(self.num_elements_total()-1)/(max_elements-1))
        alpha = np.max((min_alpha, alpha))
        if legend is False:
            legend = None
        if hasattr(self, 'history') and linewidth != 0:
            # Plot trajectories
            from matplotlib.colors import is_color_like
            if linecolor is None or is_color_like(linecolor) is True:
                if is_color_like(linecolor):
                    linecolor = linecolor
                else:
                    linecolor = 'gray'
                if compare is not None and legend is not None:
                    if legend is True:
                        if hasattr(compare, 'len'):
                            numleg = len(compare)
                        else:
                            numleg = 2
                        legend = ['Simulation %d' % (i+1) for i in
                                  range(numleg)]
                    ax.plot(x[:,0], y[:,0], color=linecolor, alpha=alpha, label=legend[0], linewidth=linewidth, transform = gcrs)
                    ax.plot(x, y, color=linecolor, alpha=alpha, label='_nolegend_', linewidth=linewidth, transform = gcrs)
                else:
                    ax.plot(x, y, color=linecolor, alpha=alpha, linewidth=linewidth, transform = gcrs)
            else:
                colorbar = True
                # Color lines according to given parameter
                try:
                    if isinstance(linecolor, str):
                        param = self.history[linecolor]
                    else:
                        param = linecolor
                except:
                    raise ValueError(
                        'Available parameters to be used for linecolors: ' +
                        str(self.history.dtype.fields))
                from matplotlib.collections import LineCollection
                for i in range(x.shape[1]):
                    vind = np.arange(index_of_first[i], index_of_last[i] + 1)
                    points = np.array(
                        [x[vind, i].T, y[vind, i].T]).T.reshape(-1, 1, 2)
                    segments = np.concatenate([points[:-1], points[1:]],
                                              axis=1)
                    if lvmin is None:
                        lvmin = param.min()
                        lvmax = param.max()
                    lc = LineCollection(segments,
                                        #cmap=plt.get_cmap('Spectral'),
                                        cmap=cmap,
                                        norm=plt.Normalize(lvmin, lvmax), transform = gcrs)
                    #lc.set_linewidth(3)
                    lc.set_array(param.T[vind, i])
                    mappable = ax.add_collection(lc)
                #axcb = fig.colorbar(lc, ax = ax, orientation = 'horizontal')
                #try:  # Add unit to colorbar if available
                #    colorbarstring = linecolor + '  [%s]' % \
                #        (self.history_metadata[linecolor]['units'])
                #except:
                #    colorbarstring = linecolor
                ##axcb.set_label(colorbarstring)
                #axcb.set_label(colorbarstring, size=14)
                #axcb.ax.tick_params(labelsize=14)

        if compare is None:
            label_initial = 'initial (%i)' % x.shape[1]
            label_active = 'active (%i)' % (x.shape[1] - self.num_elements_deactivated())
            color_initial = self.status_colors['initial']
            color_active = self.status_colors['active']
        else:
            label_initial = None
            label_active = None
            color_initial = 'gray'
            color_active = 'gray'
        if show_particles is True:
            if show_initial is True:
                ax.scatter(x[index_of_first, range(x.shape[1])],
                       y[index_of_first, range(x.shape[1])],
                       s=markersize,
                       zorder=10, edgecolor=markercolor, linewidths=.2,
                       c=color_initial, label=label_initial,
                       transform = gcrs)
            if surface_color is not None:
                color_active = surface_color
                label_active = 'surface'
            ax.scatter(x[index_of_last, range(x.shape[1])],
                       y[index_of_last, range(x.shape[1])],
                       s=markersize, zorder=3,
                       edgecolor=markercolor, linewidths=.2,
                       c=color_active, label=label_active,
                       transform = gcrs)
            #if submerged_color is not None:
            #    map.scatter(x[range(x.shape[0]), index_of_last],
            #                y[range(x.shape[0]), index_of_last], s=markersize,
            #                zorder=3, edgecolor=markercolor, linewidths=.2,
            #                c=submerged_color, label='submerged')

            x_deactivated, y_deactivated = (self.elements_deactivated.lon,
                                               self.elements_deactivated.lat)
            # Plot deactivated elements, labeled by deactivation reason
            for statusnum, status in enumerate(self.status_categories):
                if status == 'active':
                    continue  # plotted above
                if status not in self.status_colors:
                    # If no color specified, pick an unused one
                    for color in ['red', 'blue', 'green', 'black', 'gray',
                                  'cyan', 'DarkSeaGreen', 'brown']:
                        if color not in self.status_colors.values():
                            self.status_colors[status] = color
                            break
                indices = np.where(self.elements_deactivated.status == statusnum)
                if len(indices[0]) > 0:
                    if (status == 'seeded_on_land' or
                        status == 'seeded_at_nodata_position'):
                        zorder = 11
                    else:
                        zorder = 3
                    if compare is not None:
                        legstr = None
                    else:
                        legstr = '%s (%i)' % (status, len(indices[0]))
                    if compare is None:
                        color_status = self.status_colors[status]
                    else:
                        color_status = 'gray'
                    ax.scatter(x_deactivated[indices], y_deactivated[indices],
                                s=markersize,
                                zorder=zorder, edgecolor=markercolor, linewidths=.1,
                                c=color_status, label=legstr,
                                transform = gcrs)

        if compare is not None:
            cd = self._get_comparison_xy_for_plots(compare)
            for i, c in enumerate(cd):
                if legend != None:
                    legstr = legend[i+1]
                else:
                    legstr = None
                ax.plot(c['x_other'].T[:,0], c['y_other'].T[:,0], self.plot_comparison_colors[i+1] + '-', label=legstr, linewidth=linewidth, transform = gcrs)
                ax.plot(c['x_other'].T, c['y_other'].T, self.plot_comparison_colors[i+1] + '-', label='_nolegend_', linewidth=linewidth, transform = gcrs)
                ax.scatter(c['x_other'][range(c['x_other'].shape[0]), c['index_of_last_other']],
                    c['y_other'][range(c['y_other'].shape[0]), c['index_of_last_other']],
                    s=markersize,
                    zorder=3, edgecolor=markercolor, linewidths=.2,
                    c=self.plot_comparison_colors[i+1], transform = gcrs)

        try:
            if legend is not None:# and compare is None:
                plt.legend(loc=legend_loc, markerscale=2)
        except Exception as e:
            self.logger.warning('Cannot plot legend, due to bug in matplotlib:')
            self.logger.warning(traceback.format_exc())

        if background is not None:
            if hasattr(self, 'time'):
                time = self.time - self.time_step_output
            else:
                time = None
            if background == 'residence':
                scalar,lon_res,lat_res = self.get_residence_time(
                    pixelsize_m=density_pixelsize_m)
                scalar[scalar==0] = np.nan
                lon_res, lat_res = np.meshgrid(lon_res[0:-1], lat_res[0:-1])
                lon_res = lon_res.T
                lat_res = lat_res.T
                map_x, map_y = (lon_res, lat_res)
            else:
                map_x, map_y, scalar, u_component, v_component = \
                    self.get_map_background(ax, background, time=time)
                                        #self.time_step_output)

            if show_scalar is True:
                if contourlines is False:
                    scalar = np.ma.masked_invalid(scalar)
                    mappable = ax.pcolormesh(map_x, map_y, scalar, alpha=1,
                                   vmin=vmin, vmax=vmax, cmap=cmap, transform = gcrs)
                else:
                    if contourlines is True:
                        CS = ax.contour(map_x, map_y, scalar,
                                         colors='gray', transform = gcrs)
                    else:
                        # contourlines is an array of values
                        CS = ax.contour(map_x, map_y, scalar, contourlines,
                                         colors='gray', transform = gcrs)
                    plt.clabel(CS, fmt='%g')

        if mappable is not None:
            cb = fig.colorbar(mappable, orientation='horizontal', pad=.05, aspect=30, shrink=.8)
            # TODO: need better control of colorbar content
            if background is not None:
                cb.set_label(str(background))
            if linecolor is not None:
                cb.set_label(str(linecolor))

        if type(background) is list:
            delta_x = (map_x[1,2] - map_x[1,1])/2.
            delta_y = (map_y[2,1] - map_y[1,1])/2.
            ax.quiver(map_x[::skip, ::skip] + delta_x, map_y[::skip, ::skip] + delta_y,
                      u_component[::skip, ::skip],
                      v_component[::skip, ::skip], scale=scale, transform = gcrs)

        if lcs is not None:
            map_x_lcs, map_y_lcs = (lcs['lon'], lcs['lat'])
            ax.pcolormesh(
                map_x_lcs, map_y_lcs, lcs['ALCS'][0,:,:], alpha=1,
                vmin=vmin, vmax=vmax, cmap=cmap, transform = gcrs)

        if title is not None:
            if title == 'auto':
                if hasattr(self, 'time'):
                    plt.title('%s\n%s to %s UTC (%i steps)' % (
                              self._figure_title(),
                              self.start_time.strftime('%Y-%m-%d %H:%M'),
                              self.time.strftime('%Y-%m-%d %H:%M'),
                              self.steps_output))
                else:
                    plt.title('%s\n%i elements seeded at %s UTC' % (
                              self._figure_title(),
                              self.num_elements_scheduled(),
                              self.elements_scheduled_time[0].strftime(
                              '%Y-%m-%d %H:%M')))
            else:
                plt.title(title)

        if trajectory_dict is not None:
            self._plot_trajectory_dict(ax, trajectory_dict)

        #plt.gca().tick_params(labelsize=14)

        #fig.canvas.draw()
        #fig.set_tight_layout(True)
        if filename is not None:
            plt.savefig(filename)
            self.logger.info('Time to make plot: ' +
                         str(datetime.now() - start_time))
        else:
            if show is True:
                plt.show()

        return ax, plt

    def _substance_name(self):
        return None

    def _figure_title(self):
        if self._substance_name() is None:
            return 'OpenDrift - ' + type(self).__name__
        else:
            return 'OpenDrift - ' + type(self).__name__ + ' (%s)' % self._substance_name()

    def _plot_trajectory_dict(self, ax, trajectory_dict):
        '''Plot provided trajectory along with simulated'''
        time = trajectory_dict['time']
        time = np.array(time)
        i = np.where((time>=self.start_time) & (time<=self.time))[0]
        x, y = (np.atleast_1d(trajectory_dict['lon'])[i], np.atleast_1d(trajectory_dict['lat'])[i])
        ls = trajectory_dict['linestyle']

        gcrs = ccrs.PlateCarree()

        ax.plot(x, y, ls, linewidth=2, transform = gcrs)
        ax.plot(x[0], y[0], 'ok', transform = gcrs)
        ax.plot(x[-1], y[-1], 'xk', transform = gcrs)

    def get_map_background(self, ax, background, time=None):
        # Get background field for plotting on map or animation
        # TODO: this method should be made more robust
        if type(background) is list:
            variable = background[0]  # A vector is requested
        else:
            variable = background  # A scalar is requested
        for readerName in self.readers:
            reader = self.readers[readerName]
            if variable in reader.variables:
                if time is None or (time>= reader.start_time
                        and time <= reader.end_time) or (
                        reader.always_valid is True):
                    break
        if time is None:
            if hasattr(self, 'elements_scheduled_time'):
                # Using time of first seeded element
                time = self.elements_scheduled_time[0]

        # Get reader coordinates covering given map area
        axisproj = pyproj.Proj(ax.projection.proj4_params)
        xmin, xmax, ymin, ymax = ax.get_extent(ccrs.PlateCarree())
        cornerlons = np.array([xmin, xmin, xmax, xmax])
        cornerlats = np.array([ymin, ymax, ymin, ymax])
        reader_x, reader_y = reader.lonlat2xy(cornerlons, cornerlats)
        if sum(~np.isfinite(reader_x+reader_y)) > 0:
            # Axis corner points are not within reader domain
            reader_x = np.array([reader.xmin, reader.xmax])
            reader_y = np.array([reader.ymin, reader.ymax])
        else:
            reader_x = np.linspace(reader_x.min(), reader_x.max(), 10)
            reader_y = np.linspace(reader_y.min(), reader_y.max(), 10)

        data = reader.get_variables(
            background, time, reader_x, reader_y, None, block=True)
        reader_x, reader_y = np.meshgrid(data['x'], data['y'])
        if type(background) is list:
            u_component = data[background[0]]
            v_component = data[background[1]]
            with np.errstate(invalid='ignore'):
                scalar = np.sqrt(u_component**2 + v_component**2)
            # NB: rotation not completed!
            u_component, v_component = reader.rotate_vectors(
                reader_x, reader_y, u_component, v_component,
                reader.proj, ccrs.PlateCarree(
                globe=ccrs.Globe(datum='WGS84', ellipse='WGS84')
                ).proj4_init)
        else:
            scalar = data[background]
            u_component = v_component = None

        # Shift one pixel for correct plotting
        reader_x = reader_x - reader.delta_x
        reader_y = reader_y - reader.delta_y
        if reader.projected is False:
            reader_y[reader_y<0] = 0
            reader_x[reader_x<0] = 0

        rlons, rlats = reader.xy2lonlat(reader_x, reader_y)
        if rlons.max() > 360:
            rlons = rlons - 360
        map_x, map_y = (rlons, rlats)

        scalar = np.ma.masked_invalid(scalar)
        if hasattr(reader, 'convolve'):
            from scipy import ndimage
            N = reader.convolve
            if isinstance(N, (int, np.integer)):
                kernel = np.ones((N, N))
                kernel = kernel/kernel.sum()
            else:
                kernel = N
            self.logger.debug('Convolving variables with kernel: %s' % kernel)
            scalar = ndimage.convolve(
                scalar, kernel, mode='nearest')

        return map_x, map_y, scalar, u_component, v_component
    def plot_environment(self, filename=None, ax=None, show=True):
        """Plot mean wind and current velocities of element of last run."""
        x_wind = self.get_property('x_wind')[0]
        y_wind = self.get_property('y_wind')[0]
        wind = np.sqrt(x_wind**2 + y_wind**2)
        x_sea_water_velocity = self.get_property('x_sea_water_velocity')[0]
        y_sea_water_velocity = self.get_property('y_sea_water_velocity')[0]
        current = np.sqrt(x_sea_water_velocity**2 + y_sea_water_velocity**2)
        wind = np.ma.mean(wind, axis=1)
        current = np.ma.mean(current, axis=1)
        time, time_relative = self.get_time_array()
        time = np.array([t.total_seconds()/3600. for t in time_relative])

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        ax.plot(time, wind, 'b', label='Wind speed')
        ax.set_ylabel('Wind speed  [m/s]', color='b')
        ax.set_xlim([0, time[-1]])
        ax.set_ylim([0, wind.max()*1.1])

        ax2 = ax.twinx()
        ax2.plot(time, current, 'r', label='Current speed')
        ax2.set_ylabel('Current speed  [m/s]', color='r')
        ax2.set_xlim([0, time[-1]])
        ax2.set_ylim([0, current.max()*1.1])
        for tl in ax.get_yticklabels():
            tl.set_color('b')
        for tl in ax2.get_yticklabels():
            tl.set_color('r')
        ax.set_xlabel('Time  [hours]')
        ax.legend(loc='upper left')
        ax2.legend(loc='lower right')

        if filename is None:
            if show is True:
                plt.show()
        else:
            plt.savefig(filename)

    def plot_property(self, prop, mean=False):
        """Basic function to plot time series of any element properties."""
        import matplotlib.pyplot as plt
        from matplotlib import dates

        hfmt = dates.DateFormatter('%d %b %Y %H:%M')
        fig = plt.figure()
        ax = fig.gca()
        ax.xaxis.set_major_formatter(hfmt)
        plt.xticks(rotation='vertical')
        start_time = self.start_time
        # In case start_time is unsupported cftime
        start_time = datetime(start_time.year, start_time.month, start_time.day,
                              start_time.hour, start_time.minute, start_time.second)
        times = [start_time + n*self.time_step_output
                 for n in range(self.steps_output)]
        data = self.history[prop].T[0:len(times), :]
        if mean is True:  # Taking average over elements
            data = np.mean(data, axis=1)
        plt.plot(times, data)
        plt.title(prop)
        plt.xlabel('Time  [UTC]')
        try:
            plt.ylabel('%s  [%s]' %
                       (prop, self.elements.variables[prop]['units']))
        except:
            plt.ylabel(prop)
        plt.subplots_adjust(bottom=.3)
        plt.grid()
        plt.show()
    def write_geotiff(self, filename, pixelsize_km=.2):
        '''Write one GeoTiff image per timestep.

        filename should contain date identifiers, e.g. 'img_%Y%m%d_%H%M.tif'
        https://docs.python.org/2/library/datetime.html#strftime-and-strptime-behavior
        '''

        try:
            from osgeo import gdal, osr
        except:
            raise ValueError('GDAL is needed to write geotiff images.')
        import matplotlib.pyplot as plt
        driver = gdal.GetDriverByName('GTiff')
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        colortable = gdal.ColorTable()
        colortable.SetColorEntry(0,(255,255,255, 0))
        colortable.SetColorEntry(1,(0,0,0, 255))
        colortable.SetColorEntry(2,(255,0,0, 255))
        colortable.SetColorEntry(3,(0,255,0, 255))
        colortable.SetColorEntry(4,(0,0,255, 255))

        lon = self.get_property('lon')[0]
        lat = self.get_property('lat')[0]
        status = self.get_property('status')[0]
        times = self.get_time_array()[0]
        deltalat = pixelsize_km/111.0  # km to degrees
        deltalon = deltalat/np.cos(np.radians((lat.min() + lat.max())/2))
        lat_array = np.arange(lat.min()-deltalat, lat.max()+deltalat, deltalat)
        lon_array = np.arange(lon.min()-deltalat, lon.max()+deltalon, deltalon)
        ilon = (np.round((lon-lon.min())/deltalon)).astype(int)
        ilat = (np.round((lat-lat.min())/deltalat)).astype(int)
        # Setting masked values to zero, for use as indices
        ilon[ilon.mask] = 0
        ilat[ilat.mask] = 0
        status[ilon.mask] = 0
        image = np.zeros((len(times), len(lon_array),
                          len(lat_array))).astype(int)
        geotransform = [lon_array.min(), deltalon, 0,
                        lat_array.max(), 0, -deltalat]
        for i, t in enumerate(times):
            image[i, ilon[i,:], ilat[i,:]] = status[i, :] + 1
            filename_i = t.strftime(filename)
            ds = driver.Create(filename_i, len(lon_array), len(lat_array),
                               1, gdal.GDT_Byte, )
            ds.SetProjection(srs.ExportToWkt())
            ds.SetGeoTransform(geotransform)
            outband=ds.GetRasterBand(1)
            outband.SetNoDataValue(0)
            outband.WriteArray(np.fliplr(image[i, :, :]).transpose())
            outband.SetColorTable(colortable)
            ds = None



class Animator(Plotter):
    def animation(self, buffer=.2, corners=None, filename=None, compare=None,
                  background=None, bgalpha=.5, vmin=None, vmax=None, drifter=None,
                  skip=5, scale=10, color=False, clabel=None,
                  colorbar=True, cmap=None, density=False, show_elements=True,
                  show_trajectories=False, hide_landmask=False,
                  density_pixelsize_m=1000, unitfactor=1, lcs=None,
                  surface_only=False, markersize=20, origin_marker=None,
                  legend=None, legend_loc='best', fps=10, lscale=None, fast=False):
        """Animate last run."""


        if self.num_elements_total() == 0 and not hasattr(self, 'ds'):
            raise ValueError('Please run simulation before animating')

        start_time = datetime.now()
        if cmap is None:
            cmap = 'jet'
        if isinstance(cmap, str):
            cmap = matplotlib.cm.get_cmap(cmap)

        if color is False and background is None and lcs is None and density is False:
            colorbar = False

        markercolor = self.plot_comparison_colors[0]

        if isinstance(density, str):
            # Density field is weighted by this variable
            # TODO: not yet implemented!
            density_weight = density
            density = True
        else:
            if density is True:
                density_weight = None
        if density is True:  # Get density arrays
            if hasattr(self, 'ds'):  # opened with Xarray
                if origin_marker is None:
                    origin_marker = 0
                    per_origin_marker = False
                else:
                    per_origin_marker = True
                H, lon_array, lat_array = self.get_density_xarray(pixelsize_m=density_pixelsize_m,
                                            per_origin_marker=per_origin_marker)
                H = H[origin_marker]  # Presently only for origin_marker = 0
            else:
                H, H_submerged, H_stranded, lon_array, lat_array = \
                    self.get_density_array(pixelsize_m=density_pixelsize_m,
                                           weight=density_weight)
                H = H + H_submerged + H_stranded

        # x, y are longitude, latitude -> i.e. in a PlateCarree CRS
        gcrs = ccrs.PlateCarree()

        def plot_timestep(i):
            """Sub function needed for matplotlib animation."""
            ax.set_title('%s\n%s UTC' % (self._figure_title(), times[i]))
            if background is not None:
                map_x, map_y, scalar, u_component, v_component = \
                    self.get_map_background(ax, background,
                                            time=times[i])
                # https://stackoverflow.com/questions/18797175/animation-with-pcolormesh-routine-in-matplotlib-how-do-i-initialize-the-data
                bg.set_array(scalar[:-1,:-1].ravel())
                if type(background) is list:
                    bg_quiv.set_UVC(u_component[::skip, ::skip], v_component[::skip, ::skip])

            if lcs is not None:
                ax.pcolormesh(
                    lcs['lon'], lcs['lat'], lcs['ALCS'][i,:,:], alpha=bgalpha,
                    vmin=vmin, vmax=vmax, cmap=cmap, transform = gcrs)

            if density is True:
                # Update density plot
                pm.set_array(H[i,:,:].ravel())

            # Move points
            if show_elements is True:
                points.set_offsets(np.c_[x[i, range(x.shape[1])],
                                         y[i, range(x.shape[1])]])
                points_deactivated.set_offsets(np.c_[
                    x_deactive[index_of_last_deactivated < i],
                    y_deactive[index_of_last_deactivated < i]])
                if color is not False:  # Update colors
                    points.set_array(colorarray[:, i])
                    if isinstance(color, str):
                        points_deactivated.set_array(
                            colorarray_deactivated[
                                index_of_last_deactivated < i])

            if drifter is not None:
                from bisect import bisect_left
                ind = np.max(
                    (0, bisect_left(drifter['time'], times[i]) - 1))
                if i < 1 or i >= len(times)-1 or \
                        drifter['time'][ind] < times[i-1] or \
                        drifter['time'][ind] > times[i+1]:
                    # Do not show when outside time interval
                    drifter_pos.set_offsets([])
                else:
                    drifter_pos.set_offsets(
                        np.c_[drifter['x'][ind], drifter['y'][ind]])

            if show_elements is True:
                if compare is not None:
                    for cd in compare_list:
                        cd['points_other'].set_offsets(np.c_[
                            cd['x_other'][range(cd['x_other'].shape[0]), i],
                            cd['y_other'][range(cd['x_other'].shape[0]), i]])
                        cd['points_other_deactivated'].set_offsets(np.c_[
                            cd['x_other_deactive'][
                                cd['index_of_last_deactivated_other'] < i],
                            cd['y_other_deactive'][
                                cd['index_of_last_deactivated_other'] < i]])
                    return points, cd['points_other']
                else:
                    return points

        # Find map coordinates and plot points with empty data
        fig, ax, crs, x, y, index_of_first, index_of_last = \
            self.set_up_map(buffer=buffer, corners=corners, lscale=lscale,
                            fast=fast, hide_landmask=hide_landmask)

        if surface_only is True:
            z = self.get_property('z')[0]
            x[z<0] = np.nan
            y[z<0] = np.nan

        if show_trajectories is True:
            ax.plot(x, y, color='gray', alpha=.1, transform = gcrs)

        if color is not False and show_elements is True:
            if isinstance(color, str):
                colorarray = self.get_property(color)[0].T
                colorarray = colorarray*unitfactor
                colorarray_deactivated = \
                    self.get_property(color)[0][
                        index_of_last[self.elements_deactivated.ID-1],
                                      self.elements_deactivated.ID-1].T
            else:
                colorarray = color
            if vmin is None:
                vmin = colorarray.min()
                vmax = colorarray.max()

        if background is not None:
            map_x, map_y, scalar, u_component, v_component = \
                self.get_map_background(ax, background,
                                        time=self.start_time)
            bg = ax.pcolormesh(map_x, map_y, scalar[:-1,:-1], alpha=bgalpha,
                               antialiased=True, linewidth=0.0, rasterized=True,
                               vmin=vmin, vmax=vmax, cmap=cmap, transform = gcrs)
            if type(background) is list:
                bg_quiv = ax.quiver(map_x[::skip, ::skip],
                                    map_y[::skip, ::skip],
                                    u_component[::skip, ::skip],
                                    v_component[::skip, ::skip], scale=scale, transform = gcrs)

        if lcs is not None:
            if vmin is None:
                vmin = lcs['ALCS'].min()
                vmax = lcs['ALCS'].max()
            lcsh = ax.pcolormesh(lcs['lon'], lcs['lat'], lcs['ALCS'][0,:,:],
                                  vmin=vmin, vmax=vmax, cmap=cmap,
                                  transform = gcrs)

        times = self.get_time_array()[0]
        if show_elements is True:
            index_of_last_deactivated = \
                index_of_last[self.elements_deactivated.ID-1]
        if legend is None:
            legend = ['']

        if color is False:
            c = markercolor
        else:
            c = []
        points = ax.scatter([], [], c=c, zorder=10,
                            edgecolor=[], cmap=cmap, s=markersize,
                            vmin=vmin, vmax=vmax, label=legend[0], transform = gcrs)
        # Plot deactivated elements, with transparency
        points_deactivated = ax.scatter([], [], c=c, zorder=9,
                                        vmin=vmin, vmax=vmax, s=markersize, cmap=cmap,
                                        edgecolor=[], alpha=.3, transform = gcrs)
        x_deactive, y_deactive = (self.elements_deactivated.lon, self.elements_deactivated.lat)

        if compare is not None:
            compare_list = self._get_comparison_xy_for_plots(compare)

            for cn, cd in enumerate(compare_list):
                if legend != ['']:
                    legstr = legend[cn+1]
                else:
                    legstr = None
                cd['points_other'] = \
                    ax.scatter([], [], c=self.plot_comparison_colors[cn+1],
                               s=markersize, label=legstr, zorder=10, transform = gcrs)
                # Plot deactivated elements, with transparency
                cd['points_other_deactivated'] = \
                    ax.scatter([], [], alpha=.3, zorder=9,
                               c=self.plot_comparison_colors[cn+1],
                               s=markersize, transform = gcrs)

            if legend != ['', '']:
                plt.legend(markerscale=2, loc=legend_loc)

        if density is True:
            cmap.set_under('w')
            H = np.ma.masked_where(H==0, H)
            lat_array, lon_array = np.meshgrid(lat_array, lon_array)
            if vmax is None:
                vmax=H.max()
            pm = ax.pcolormesh(lon_array, lat_array, H[0,:,:],
                               vmin=0.1, vmax=vmax, cmap=cmap, transform=gcrs)

        if drifter is not None:
            drifter['x'], drifter['y'] = (drifter['lon'], drifter['lat'])
            #map.plot(drifter['x'], drifter['y'])
            drifter_pos = ax.scatter([], [], c='r', zorder=15,
                                     label='Drifter', transform = gcrs)

        fig.canvas.draw()
        fig.set_tight_layout(True)
        if colorbar is True:
            if color is not False:
                if isinstance(color, str) or clabel is not None:
                    if clabel is None:
                        clabel = color
                item = points
            elif density is not False:
                item = pm
                if clabel is None:
                    clabel = 'density'
            elif lcs is not None:
                item = lcsh
                if clabel is None:
                    clabel = 'LCS'
            elif background is not None:
                item = bg
                if clabel is None:
                    clabel = background

            cb = fig.colorbar(item, orientation='horizontal', pad=.05, aspect=30, shrink=.8)
            cb.set_label(clabel)
            cb.set_alpha(1)
            cb.draw_all()

        anim = animation.FuncAnimation(
            plt.gcf(), plot_timestep, blit=False,
            frames=x.shape[0], interval=50)


        if filename is not None or 'sphinx_gallery' in sys.modules:
            self._save_animation(anim, filename, fps)
            self.logger.debug('Time to make animation: %s' %
                          (datetime.now()-start_time))
        else:
            try:
                plt.show()
            except AttributeError:
                pass

    def animation_profile(self, filename=None, compare=None,
                          legend=['', ''], markersize=5, fps=20):
        """Animate vertical profile of the last run."""

        def plot_timestep(i):
            """Sub function needed for matplotlib animation."""
            #plt.gcf().gca().set_title(str(i))
            ax.set_title('%s UTC' % times[i])
            points.set_data(x[range(x.shape[0]), i],
                            z[range(x.shape[0]), i])
            points_deactivated.set_data(
                x_deactive[index_of_last_deactivated < i],
                z_deactive[index_of_last_deactivated < i])

            if compare is not None:
                points_other.set_data(x_other[range(x_other.shape[0]), i],
                                      z_other[range(x_other.shape[0]), i])
                points_other_deactivated.set_data(
                    x_other_deactive[index_of_last_deactivated_other < i],
                    z_other_deactive[index_of_last_deactivated_other < i])
                return points, points_other
            else:
                return points

        # Set up plot
        index_of_first, index_of_last = \
            self.index_of_activation_and_deactivation()
        z = self.get_property('z')[0].T
        x = self.get_property('lon')[0].T
        #seafloor_depth = \
        #    -self.get_property('sea_floor_depth_below_sea_level')[0].T
        fig = plt.figure(figsize=(10, 6.))  # Suitable aspect ratio
        ax = fig.gca()
        plt.xlabel('Longitude [degrees]')
        plt.ylabel('Depth [m]')
        times = self.get_time_array()[0]
        index_of_last_deactivated = \
            index_of_last[self.elements_deactivated.ID-1]
        points = plt.plot([], [], '.k', label=legend[0],
                          markersize=markersize)[0]
        # Plot deactivated elements, with transparency
        points_deactivated = plt.plot([], [], '.k', alpha=.3)[0]
        x_deactive = self.elements_deactivated.lon
        z_deactive = self.elements_deactivated.z

        if compare is not None:
            if type(compare) is str:
                # Other is given as filename
                other = self.__class__(loglevel=0)
                other.io_import_file(compare)
            else:
                # Other is given as an OpenDrift object
                other = compare
            z_other = other.get_property('z')[0].T
            x_other = self.get_property('lon')[0].T
            points_other = plt.plot(x_other[0, 0], z_other[0, 0], '.r',
                                    label=legend[1],
                                    markersize=markersize)[0]
            # Plot deactivated elements, with transparency
            points_other_deactivated = plt.plot([], [], '.r', alpha=.3)[0]
            x_other_deactive = other.elements_deactivated.lon
            z_other_deactive = other.elements_deactivated.z
            firstlast = np.ma.notmasked_edges(x_other, axis=1)
            index_of_last_other = firstlast[1][1]
            index_of_last_deactivated_other = \
                index_of_last_other[other.elements_deactivated.ID-1]
            xmax = np.maximum(x.max(), x_other.max())
            xmin = np.minimum(x.min(), x_other.min())
            zmax = np.maximum(z.max(), z_other.max())
            zmin = np.minimum(z.min(), z_other.min())
        else:
            xmin = x.min()
            xmax = x.max()
            zmin = z.min()
            zmax = z.max()

        # Set figure limits
        sky = (zmax-zmin)*.1  # Sky height is 10% of water depth
        plt.xlim([xmin, xmax])
        plt.ylim([zmin, sky])
        ax.add_patch(plt.Rectangle((xmin, 0), xmax-xmin, sky,
                     color='lightsteelblue'))
        ax.add_patch(plt.Rectangle((xmin, zmin), xmax-xmin,
                     -zmin, color='cornflowerblue'))

        if legend != ['', '']:
            plt.legend(loc=4)

        anim = animation.FuncAnimation(plt.gcf(), plot_timestep, blit=False,
                                       frames=x.shape[1], interval=150)

        if filename is not None or 'sphinx_gallery' in sys.modules:
            self._save_animation(anim, filename, fps)

        else:
            try:
                plt.show()
            except AttributeError:
                pass

    def _save_animation(self, anim, filename, fps):
        if 'sphinx_gallery' in sys.modules:
            # This assumes that the calling script is two frames up in the stack. If
            # _save_animation is called through a more deeply nested method, it will
            # not give the correct result.
            caller = inspect.stack()[2]
            caller = os.path.splitext(os.path.basename(caller.filename))[0]

            # Calling script is string input (e.g. from ..plot::)
            if caller == '<string>':
                caller = 'plot_directive'
                adir = os.path.realpath('../source/gallery/animations')
            else:
                adir = os.path.realpath('../docs/source/gallery/animations')

            if not hasattr(OpenDriftSimulation, '__anim_no__'):
                OpenDriftSimulation.__anim_no__ = { }

            if caller not in  OpenDriftSimulation.__anim_no__:
                OpenDriftSimulation.__anim_no__[caller] = 0

            os.makedirs(adir, exist_ok=True)

            filename = '%s_%d.gif' % (caller, OpenDriftSimulation.__anim_no__[caller])
            OpenDriftSimulation.__anim_no__[caller] += 1

            filename = os.path.join(adir, filename)

        self.logger.info('Saving animation to ' + filename + '...')

        try:
            if filename[-4:] == '.gif':  # GIF
                self.logger.info('Making animated gif...')
                anim.save(filename, fps=fps, writer='imagemagick')
            else:  # MP4
                try:
                    try:
                        # For perfect quality, but larger file size
                        #FFwriter=animation.FFMpegWriter(fps=fps, extra_args=['-vcodec', 'libx264'])
                        FFwriter=animation.FFMpegWriter(fps=fps,
                            codec='libx264', bitrate=1800,
                            extra_args=['-profile:v', 'baseline',
                                        '-pix_fmt', 'yuv420p', '-an'])
                        anim.save(filename, writer=FFwriter)
                    except Exception as e:
                        self.logger.info(e)
                        anim.save(filename, fps=fps, bitrate=1800,
                                  extra_args=['-an',
                                              '-pix_fmt', 'yuv420p'])
                except Exception as e:
                    self.logger.info(e)
                    anim.save(filename, fps=fps)
                    self.logger.warning('Animation might not be HTML5 compatible.')

        except Exception as e:
            self.logger.info('Could not save animation:')
            self.logger.info(e)
            self.logger.debug(traceback.format_exc())

        if 'sphinx_gallery' in sys.modules:
            plt.close()

