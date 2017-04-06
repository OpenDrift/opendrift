#!/usr/bin/env python

import sys
from datetime import datetime, timedelta
import numpy as np

from PIL import ImageTk, Image
import Tkinter as tk
from opendrift.models.openoil3D import OpenOil3D
from opendrift.models.leeway import Leeway
from opendrift.readers import reader_netCDF_CF_generic
from opendrift.readers import reader_basemap_landmask


class TextRedirector(object):
    def __init__(self, widget, tag="stdout"):
        self.defstdout = sys.stdout
        self.widget = widget
        self.tag = tag

    def write(self, str):
        self.widget.configure(state="normal")
        self.widget.insert("end", str, (self.tag,))
        self.widget.update_idletasks()
        self.widget.see(tk.END)

    def flush(self):
        self.defstdout.flush()


class OpenDriftGUI(tk.Tk):

    def __init__(self):
        tk.Tk.__init__(self)
        self.title('OpenDrift')
        o = OpenOil3D(weathering_model='noaa', location='NORWAY')
        try:
            img = ImageTk.PhotoImage(Image.open(o.test_data_folder() +
                                     '../../docs/opendrift_logo.png'))
            panel = tk.Label(self.master, image=img)
            panel.image = img
            panel.grid(row=0, column=0)
        except:
            pass # Could not display logo

        self.top = tk.Frame(self.master,
                            relief=tk.FLAT, pady=25, padx=25)
        self.top.grid(row=0, column=1, rowspan=1)

        tk.Label(self.top, text='Simulation type').grid(row=0, column=0)
        self.model = tk.StringVar()
        models = ['OpenOil', 'Leeway']
        self.model.set(models[0])
        self.modeldrop = tk.OptionMenu(self.top, self.model,
                                       *(models), command=self.set_model)
        self.modeldrop.grid(row=0, column=1)

        self.categoryLabel = tk.Label(self.master, text='Oil type')
        self.categoryLabel.grid(row=1, column=0)
        oljetyper = o.oiltypes
        self.oljetype = tk.StringVar()
        self.oljetype.set(oljetyper[0])
        self.categorydrop = tk.OptionMenu(self.master,
                                          self.oljetype, *oljetyper)
        self.categorydrop.grid(row=1, column=1)

        ##########
        # Release
        ##########
        self.start_t = tk.Frame(self.master, relief=tk.FLAT)
        self.start_t.grid(row=2, column=0, rowspan=1)
        startlabel = tk.Label(self.start_t, text="\n\nStart release\n\n")
        startlabel.grid(row=0, column=0)
        self.start = tk.Frame(self.master, bg='lightgray', bd=2,
                              relief=tk.SUNKEN, pady=5, padx=5)
        self.start.grid(row=2, column=1, rowspan=1)

        tk.Label(self.start, text='Longitude').grid(row=0, column=0)
        tk.Label(self.start, text='Latitude').grid(row=0, column=1)
        tk.Label(self.start, text='Radius [m]').grid(row=0, column=2)
        self.lon = tk.Entry(self.start, width=6, justify=tk.RIGHT)
        self.lat = tk.Entry(self.start, width=6, justify=tk.RIGHT)
        self.radius = tk.Entry(self.start, width=6, justify=tk.RIGHT)
        self.lon.grid(row=1, column=0)
        self.lon.insert(0, '4.5')
        self.lat.grid(row=1, column=1)
        self.lat.insert(0, '60.0')
        self.radius.grid(row=1, column=2)
        self.radius.insert(0, '5000')
        ##########
        # Time
        ##########
        now = datetime.now()
        tk.Label(self.start, text='Day').grid(row=2, column=0)
        tk.Label(self.start, text='Month').grid(row=2, column=1)
        tk.Label(self.start, text='Year').grid(row=2, column=2)
        tk.Label(self.start, text='Hour').grid(row=2, column=3)
        tk.Label(self.start, text='Minutes [UTC]').grid(row=2, column=4)
        self.datevar = tk.StringVar()
        self.dates = range(1, 32)
        self.datevar.set(now.day)
        self.date = tk.OptionMenu(self.start, self.datevar, *self.dates)
        self.date.grid(row=3, column=0)

        self.monthvar = tk.StringVar()
        self.months = ['January', 'February', 'March', 'April', 'May',
                       'June', 'July', 'August', 'September', 'October',
                       'November', 'December']
        self.monthvar.set(self.months[now.month-1])
        self.month = tk.OptionMenu(self.start, self.monthvar,
                                   *self.months)
        self.month.grid(row=3, column=1)

        self.yearvar = tk.StringVar()
        self.years = range(2015, 2017)
        self.yearvar.set(now.year)
        self.year = tk.OptionMenu(self.start, self.yearvar, *self.years)
        self.year.grid(row=3, column=2)

        self.hourvar = tk.StringVar()
        self.hours = range(0, 24)
        self.hourvar.set(now.hour)
        self.hour = tk.OptionMenu(self.start, self.hourvar, *self.hours)
        self.hour.grid(row=3, column=3)

        self.minutevar = tk.StringVar()
        self.minutes = range(0, 60, 5)
        self.minutevar.set(now.minute)
        self.minute = tk.OptionMenu(self.start, self.minutevar,
                                    *self.minutes)
        self.minute.grid(row=3, column=4)

        ###############
        # Release End
        ###############
        self.end_t = tk.Frame(self.master, relief=tk.FLAT)
        endlabel = tk.Label(self.end_t, text="\n\nEnd release\n\n")
        endlabel.grid(row=0, column=0)
        self.end_t.grid(row=3, column=0, rowspan=1)

        self.end = tk.Frame(self.master, bg='gray', bd=2,
                            relief=tk.SUNKEN, padx=5, pady=5)
        self.end.grid(row=3, column=1)

        tk.Label(self.end, text='Longitude').grid(row=0, column=0)
        tk.Label(self.end, text='Latitude').grid(row=0, column=1)
        tk.Label(self.end, text='Radius [m]').grid(row=0, column=2)
        self.elon = tk.Entry(self.end, width=6, justify=tk.RIGHT)
        self.elat = tk.Entry(self.end, width=6, justify=tk.RIGHT)
        self.eradius = tk.Entry(self.end, width=6, justify=tk.RIGHT)
        self.elon.grid(row=1, column=0)
        self.elon.insert(0, '4.5')
        self.elat.grid(row=1, column=1)
        self.elat.insert(0, '60.0')
        self.eradius.grid(row=1, column=2)
        self.eradius.insert(0, '0')
        ##########
        # Time
        ##########
        now = datetime.now()
        tk.Label(self.end, text='Day').grid(row=2, column=0)
        tk.Label(self.end, text='Month').grid(row=2, column=1)
        tk.Label(self.end, text='Year').grid(row=2, column=2)
        tk.Label(self.end, text='Hour').grid(row=2, column=3)
        tk.Label(self.end, text='Minutes [UTC]').grid(row=2, column=4)
        self.edatevar = tk.StringVar()
        self.edates = range(1, 32)
        self.edatevar.set(now.day)
        self.edate = tk.OptionMenu(self.end, self.edatevar, *self.edates)
        self.edate.grid(row=3, column=0)

        self.emonthvar = tk.StringVar()
        self.emonthvar.set(self.months[now.month-1])
        self.emonth = tk.OptionMenu(self.end, self.emonthvar,
                                    *self.months)
        self.emonth.grid(row=3, column=1)

        self.eyearvar = tk.StringVar()
        self.eyears = range(2015, 2017)
        self.eyearvar.set(now.year)
        self.eyear = tk.OptionMenu(self.end, self.eyearvar, *self.eyears)
        self.eyear.grid(row=3, column=2)

        self.ehourvar = tk.StringVar()
        self.ehours = range(0, 24)
        self.ehourvar.set(now.hour)
        self.ehour = tk.OptionMenu(self.end, self.ehourvar, *self.ehours)
        self.ehour.grid(row=3, column=3)

        self.eminutevar = tk.StringVar()
        self.eminutes = range(0, 60, 5)
        self.eminutevar.set(now.minute)
        self.eminute = tk.OptionMenu(self.end, self.eminutevar,
                                     *self.eminutes)
        self.eminute.grid(row=3, column=4)

        #######################
        # Simulation duration
        #######################
        self.duration = tk.Frame(self.master, bd=2,
                                 relief=tk.FLAT, padx=5, pady=15)
        self.duration.grid(row=4, column=1)
        tk.Label(self.duration, text='Run simulation ').grid(row=4, column=0)
        self.durationhours = tk.Entry(self.duration, width=3,
                                      justify=tk.RIGHT)
        self.durationhours.grid(row=4, column=1)
        self.durationhours.insert(0, 12)
        tk.Label(self.duration, text=' hours ').grid(row=4, column=2)

        self.directionvar = tk.StringVar()
        self.directionvar.set('forwards')
        self.direction = tk.OptionMenu(self.duration, self.directionvar,
                                       'forwards', 'backwards')
        self.direction.grid(row=4, column=3)
        tk.Label(self.duration, text=' in time ').grid(row=4, column=4)

        ##############
        # Output box
        ##############
        self.text = tk.Text(self, wrap="word", height=18)
        self.text.grid(row=5, columnspan=5, sticky='nsew')
        self.text.tag_configure("stderr", foreground="#b22222")
        sys.stdout = TextRedirector(self.text, "stdout")
        sys.stderr = TextRedirector(self.text, "stderr")
        s = tk.Scrollbar(self)
        s.grid(row=5, column=5, sticky='ns')
        s.config(command=self.text.yview)
        self.text.config(yscrollcommand=s.set)

        ##############
        # Driver data
        ##############
        o = OpenOil3D()
        self.current = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be')
        #    o.test_data_folder() +
        #    '16Nov2015_NorKyst_z_surface/norkyst800_subset_16Nov2015.nc')
        self.wind = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')
        #self.wind = reader_netCDF_CF_generic.Reader(o.test_data_folder() +
        #    '16Nov2015_NorKyst_z_surface/arome_subset_16Nov2015.nc')

        print '#'*41
        print 'Current data coverage:'
        print str(self.current.start_time) + ' - ' + \
            str(self.current.end_time)
        print '#'*41
        print 'Wind data coverage:'
        print str(self.wind.start_time) + ' - ' + \
            str(self.wind.end_time)
        print '#'*41
        self.start_time = self.current.start_time

        ##########
        # RUN
        ##########
        tk.Button(self.master, text='PEIS PAO', bg='green',
                  command=self.run_opendrift).grid(row=7, column=1,
                                              sticky=tk.W, pady=4)

    def set_model(self, model):
        if model == 'OpenOil':
            self.categoryLabel['text'] = 'Oil type'
            self.oljetype.set('')
            self.categorydrop['menu'].delete(0, 'end')
            self.o = OpenOil3D(weathering_model='noaa', location='NORWAY')
            for cat in self.o.oiltypes:
                self.categorydrop['menu'].add_command(
                label=cat, command=tk._setit(self.oljetype, cat))
            self.oljetype.set(self.o.oiltypes[0])

        if model == 'Leeway':
            self.categoryLabel['text'] = 'Object type'
            self.oljetype.set('')
            self.categorydrop['menu'].delete(0, 'end')
            self.o = Leeway()
            self.leewaycategories = [
                self.o.leewayprop[c]['Description'].strip().
                replace('>', '') for c in self.o.leewayprop]
            for cat in self.leewaycategories:
                self.categorydrop['menu'].add_command(
                label=cat, command=tk._setit(self.oljetype, cat))
            self.oljetype.set(self.leewaycategories[0])

    def run_opendrift(self):
        sys.stdout.write('running OpenDrift')
        month = np.int(self.months.index(self.monthvar.get()) + 1)
        start_time = datetime(np.int(self.yearvar.get()), month,
                              np.int(self.datevar.get()),
                              np.int(self.hourvar.get()),
                              np.int(self.minutevar.get()))
        if start_time > self.current.end_time:
            sys.stdout.write('Start time after end of current data!')
            start_time = self.current.start_time
        emonth = np.int(self.months.index(self.emonthvar.get()) + 1)
        end_time = datetime(np.int(self.eyearvar.get()), emonth,
                            np.int(self.edatevar.get()),
                            np.int(self.ehourvar.get()),
                            np.int(self.eminutevar.get()))
        sys.stdout.flush()
        lon = np.float(self.lon.get())
        lat = np.float(self.lat.get())
        radius = np.float(self.radius.get())
        elon = np.float(self.elon.get())
        elat = np.float(self.elat.get())
        eradius = np.float(self.eradius.get())
        if lon != elon or lat != elat or start_time != end_time:
            lon = [lon, elon]
            lat = [lat, elat]
            radius = [radius, eradius]
            start_time = [start_time, end_time]
            cone = True
        else:
            cone = False
        if self.model.get() == 'Leeway':
            o = Leeway(loglevel=0)
            for ln, lc in enumerate(self.leewaycategories):
                if self.oljetype.get() == lc.strip().replace('>', ''):
                    print 'Leeway object category: ' + lc
                    break
            o.seed_elements(lon=lon, lat=lat, number=2000,
                            radius=radius, time=start_time,
                            objectType=ln + 1)
        if self.model.get() == 'OpenOil':
            o = OpenOil3D(weathering_model='noaa', loglevel=0)
            o.seed_elements(lon=lon, lat=lat, number=2000, radius=radius,
                            time=start_time, cone=cone,
                            oiltype=self.oljetype.get())

        o.add_reader([self.current, self.wind])
        o.set_config('general:basemap_resolution', 'h')

        time_step = 1800  # Half hour
        duration = int(self.durationhours.get())*3600/time_step
        if self.directionvar.get() == 'backwards':
            time_step = -time_step
        o.run(steps=duration, time_step=time_step)
        print o

        diana_filename = 'diana.nc'
        tk.Button(self.master, text='Save to Diana',
                  command=o.write_netcdf_density_map(diana_filename)
                  ).grid(row=7, column=2, sticky=tk.W, pady=4)
        tk.Button(self.master, text='Animation',
                  command=o.animation).grid(row=7, column=3,
                                            sticky=tk.W, pady=4)
        if self.model.get() == 'OpenOil':
            tk.Button(self.master, text='Oil Budget',
                      command=o.plot_oil_budget).grid(row=7, column=4,
                                                      sticky=tk.W, pady=4)

        #o.plot()


if __name__ == '__main__':
    OpenDriftGUI().mainloop()
