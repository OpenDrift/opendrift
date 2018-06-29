#!/usr/bin/env python

import sys
import os
from datetime import datetime, timedelta
import numpy as np

from PIL import ImageTk, Image
import Tkinter as tk
import ttk
import opendrift.models
from opendrift.models.openoil3D import OpenOil3D
from opendrift.models.leeway import Leeway
from opendrift.models.shipdrift import ShipDrift
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

        # Available models
        self.available_models = opendrift.get_model_names()

        tk.Tk.__init__(self)
        
        ##################
        # Layout frames
        ##################
        self.n = ttk.Notebook(self.master)
        self.n.grid()
        self.seed = ttk.Frame(self.n)
        self.config = ttk.Frame(self.n)
        self.n.add(self.seed, text='Seeding')
        self.n.add(self.config, text='Config')

        # Top
        self.top = tk.Frame(self.seed,
                            relief=tk.FLAT, pady=25, padx=25)
        self.top.grid(row=0, column=1, rowspan=1)
        # Config
        self.con = tk.Label(self.config, text="\n\nConfiguration\n\n")
        self.con.grid(row=0, column=1, rowspan=1)
        # Time start and end
        self.start_t = tk.Frame(self.seed, relief=tk.FLAT)
        self.start_t.grid(row=2, column=0, rowspan=1)
        self.end_t = tk.Frame(self.seed, relief=tk.FLAT)
        self.end_t.grid(row=3, column=0, rowspan=1)
        self.start = tk.Frame(self.seed, bg='lightgray', bd=2,
                              relief=tk.SUNKEN, pady=5, padx=5)
        self.start.grid(row=2, column=1, rowspan=1)
        self.end = tk.Frame(self.seed, bg='gray', bd=2,
                            relief=tk.SUNKEN, padx=5, pady=5)
        self.end.grid(row=3, column=1)
        self.coastline = tk.Frame(self.seed, bd=2,
                                 relief=tk.FLAT, padx=5, pady=0)
        self.coastline.grid(row=4, column=1)
        self.duration = tk.Frame(self.seed, bd=2,
                                 relief=tk.FLAT, padx=5, pady=5)
        self.duration.grid(row=5, column=1)
        self.seed_frame = tk.Frame(self.seed, bd=2,
                                   relief=tk.FLAT, padx=5, pady=0)
        self.seed_frame.grid(row=4, columnspan=8, sticky='nsew')
        self.output = tk.Frame(self.seed, bd=2,
                               relief=tk.FLAT, padx=5, pady=0)
        self.output.grid(row=7, column=2, columnspan=8, sticky='nsew')


        ##########################
        self.title('OpenDrift')
        o = OpenOil3D(weathering_model='noaa', location='NORWAY')
        try:
            img = ImageTk.PhotoImage(Image.open(o.test_data_folder() +
                                     '../../docs/opendrift_logo.png'))
            panel = tk.Label(self.seed, image=img)
            panel.image = img
            panel.grid(row=0, column=0)
        except:
            pass # Could not display logo
        #######################################################
        tk.Label(self.top, text='Simulation type').grid(row=0, column=0)
        self.model = tk.StringVar()
        models = ['OpenOil', 'Leeway']
        models = opendrift.get_model_names()
        self.model.set(models[0])
        self.modeldrop = tk.OptionMenu(self.top, self.model,
                                       *(models), command=self.set_model)
        self.modeldrop.grid(row=0, column=1)

        help_button = tk.Button(self.top, text='Help',
                                command=self.show_help)
        help_button.grid(row=0, column=2, padx=50)

        self.categoryLabel = tk.Label(self.seed, text='Oil type')
        self.categoryLabel.grid(row=1, column=0)
        oljetyper = o.oiltypes
        self.oljetype = tk.StringVar()
        self.oljetype.set(oljetyper[0])
        self.categorydrop = ttk.Combobox(self.seed, width=50,
                                         textvariable=self.oljetype,
                                         values=oljetyper)
        self.categorydrop.grid(row=1, column=1)

        ##########
        # Release
        ##########
        startlabel = tk.Label(self.start_t, text="\n\nStart release\n\n")
        startlabel.grid(row=0, column=0)

        tk.Label(self.start, text='Longitude').grid(row=0, column=1)
        tk.Label(self.start, text='Latitude').grid(row=0, column=0)
        tk.Label(self.start, text='Radius [m]').grid(row=0, column=2)
        self.latvar = tk.StringVar()
        self.lonvar = tk.StringVar()
        self.radiusvar = tk.StringVar()
        self.lat = tk.Entry(self.start, textvariable=self.latvar,
                            width=6, justify=tk.RIGHT)
        self.lon = tk.Entry(self.start, textvariable=self.lonvar,
                            width=6, justify=tk.RIGHT)
        self.radius = tk.Entry(self.start, width=6,
                               textvariable=self.radiusvar,
                               justify=tk.RIGHT)
        self.lon.grid(row=1, column=1)
        self.lon.insert(0, '4.5')
        self.lat.grid(row=1, column=0)
        self.lat.insert(0, '60.0')
        self.radius.grid(row=1, column=2)
        self.radius.insert(0, '1000')
        self.lonvar.trace('w', self.copy_position)
        self.latvar.trace('w', self.copy_position)
        self.radiusvar.trace('w', self.copy_position)
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
        self.years = range(2015, now.year+1)
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

        self.datevar.trace('w', self.copy_position)
        self.monthvar.trace('w', self.copy_position)
        self.yearvar.trace('w', self.copy_position)
        self.hourvar.trace('w', self.copy_position)
        self.minutevar.trace('w', self.copy_position)

        ###############
        # Release End
        ###############
        endlabel = tk.Label(self.end_t, text="\n\nEnd release\n\n")
        endlabel.grid(row=0, column=0)
        tk.Label(self.end, text='Longitude', bg='gray').grid(row=0, column=1)
        tk.Label(self.end, text='Latitude', bg='gray').grid(row=0, column=0)
        tk.Label(self.end, text='Radius [m]', bg='gray').grid(row=0, column=2)
        self.elat = tk.Entry(self.end, width=6, justify=tk.RIGHT)
        self.elon = tk.Entry(self.end, width=6, justify=tk.RIGHT)
        self.eradius = tk.Entry(self.end, width=6, justify=tk.RIGHT)
        self.elon.grid(row=1, column=1)
        self.elon.insert(0, '4.5')
        self.elat.grid(row=1, column=0)
        self.elat.insert(0, '60.0')
        self.eradius.grid(row=1, column=2)
        self.eradius.insert(0, '1000')
        ##########
        # Time
        ##########
        now = datetime.now()
        tk.Label(self.end, text='Day', bg='gray').grid(row=2, column=0)
        tk.Label(self.end, text='Month', bg='gray').grid(row=2, column=1)
        tk.Label(self.end, text='Year', bg='gray').grid(row=2, column=2)
        tk.Label(self.end, text='Hour', bg='gray').grid(row=2, column=3)
        tk.Label(self.end, text='Minutes [UTC]', bg='gray').grid(row=2, column=4)
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
        self.eyears = range(2015, now.year+1)
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
        self.eyear.config(bg='gray')
        self.emonth.config(bg='gray')
        self.edate.config(bg='gray')
        self.ehour.config(bg='gray')
        self.eminute.config(bg='gray')

        # Check seeding
        check_seed = tk.Button(self.end_t, text='Check seeding',
                                command=self.check_seeding)
        check_seed.grid(row=1, column=0, padx=0)

        #######################
        # Simulation duration
        #######################
        tk.Label(self.coastline, text='Coastline resolution ').grid(
                 row=4, column=1)
        self.mapresvar = tk.StringVar()
        self.mapres = tk.OptionMenu(self.coastline, self.mapresvar,
                                    *['full', 'high'])
        self.mapres.grid(row=4, column=2)
        self.mapresvar.set('high')
        
        tk.Label(self.duration, text='Run simulation ').grid(row=5, column=0)
        self.durationhours = tk.Entry(self.duration, width=3,
                                      justify=tk.RIGHT)
        self.durationhours.grid(row=5, column=1)
        self.durationhours.insert(0, 12)
        tk.Label(self.duration, text=' hours ').grid(row=5, column=2)

        self.directionvar = tk.StringVar()
        self.directionvar.set('forwards')
        self.direction = tk.OptionMenu(self.duration, self.directionvar,
                                       'forwards', 'backwards')
        self.direction.grid(row=5, column=3)
        tk.Label(self.duration, text=' in time ').grid(row=5, column=4)

        ##############
        # Output box
        ##############
        self.text = tk.Text(self.output, wrap="word", height=18)
        self.text.grid(row=6, columnspan=8, sticky='nsew')
        self.text.tag_configure("stderr", foreground="#b22222")
        sys.stdout = TextRedirector(self.text, "stdout")
        sys.stderr = TextRedirector(self.text, "stderr")
        s = tk.Scrollbar(self)
        s.grid(row=6, column=8, sticky='ns')
        s.config(command=self.text.yview)
        self.text.config(yscrollcommand=s.set)

        # Diana
        self.dianadir = '/vol/vvfelles/opendrift/output/'
        if os.path.exists(self.dianadir):
            self.has_diana = True
            print('Diana is available!')
            self.outputdir = '/vol/vvfelles/opendrift/output_native/'
        else:
            self.has_diana = False

        ##############
        # Initialise
        ##############
        o = OpenOil3D()

        for m in self.available_models:
            print(m)
            print(opendrift.get_model(m))
            print(opendrift.get_model(m)())

        ##########
        # RUN
        ##########
        tk.Button(self.seed, text='PEIS PAO', bg='green',
                  command=self.run_opendrift).grid(row=8, column=1,
                                              sticky=tk.W, pady=4)

    def copy_position(self, a, b, c):
        self.elat.delete(0, tk.END)
        self.elat.insert(0, self.lat.get())
        self.elon.delete(0, tk.END)
        self.elon.insert(0, self.lon.get())
        self.eradius.delete(0, tk.END)
        self.eradius.insert(0, self.radius.get())
        self.edatevar.set(self.datevar.get())
        self.emonthvar.set(self.monthvar.get())
        self.eyearvar.set(self.yearvar.get())
        self.ehourvar.set(self.hourvar.get())
        self.eminutevar.set(self.minutevar.get())

    def set_model(self, model):
        if model == 'OpenOil':
            self.categoryLabel['text'] = 'Oil type'
            self.oljetype.set('')
            self.o = OpenOil3D(weathering_model='noaa', location='NORWAY')
            self.categorydrop['values'] = self.o.oiltypes
            self.oljetype.set(self.o.oiltypes[0])
        elif model == 'Leeway':
            self.categoryLabel['text'] = 'Object type'
            self.oljetype.set('')
            self.o = Leeway()
            self.leewaycategories = [
                self.o.leewayprop[c]['Description'].strip().
                replace('>', '') for c in self.o.leewayprop]
            self.categorydrop['values'] = self.leewaycategories
            self.oljetype.set(self.leewaycategories[0])
        elif model == 'ShipDrift':
            self.categoryLabel['text'] = 'Object type'
            self.oljetype.set('')
            self.o = ShipDrift()

        print(dir(self.config))
        for con in self.config.winfo_children():
            con.destroy()
        self.con = tk.Label(self.config, text="\n\nConfiguration\n\n")
        self.con.grid(row=0, column=1, rowspan=1)
        print(dir(self.config))
        for i, cs in enumerate(self.o._config_hashstrings()):
            tk.Label(self.config, text=cs).grid(row=i, column=1, rowspan=1)

        print('Setting model: ' + model)
        print(self.o.list_configspec())
        #default, minimum, maximum, options = self.o.get_configspec_default_min_max('seed:oil_type')
        #print default
        sc = self.o.get_seed_config()
        print(sc)
        self.seed_input = {}
        self.seed_input_var = {}
        self.seed_input_label = {}
        self.seed_frame.destroy()
        self.seed_frame = tk.Frame(self.seed, bd=2,
                                   relief=tk.FLAT, padx=5, pady=0)
        self.seed_frame.grid(row=4, columnspan=8, sticky='nsew')
        # FIND
        for num, i in enumerate(sc):
            self.seed_input_label[i] = tk.Label(self.seed_frame,
                                                text=i + '\t')
            self.seed_input_label[i].grid(row=num, column=0)
            self.seed_input_var[i] = tk.StringVar()
            if type(sc[i]['options']) is list:
                self.seed_input[i] = ttk.Combobox(
                    self.seed_frame, width=50,
                    textvariable=self.seed_input_var[i],
                    values=sc[i]['options'])
                self.seed_input_var[i].set(sc[i]['default'])
            else:
                self.seed_input[i] = tk.Entry(
                    self.seed_frame, textvariable=self.seed_input_var[i],
                    width=6, justify=tk.RIGHT)
                self.seed_input[i].insert(0, sc[i]['default'])
            self.seed_input[i].grid(row=num, column=1)

    def show_help(self):
        help_url = 'https://github.com/OpenDrift/opendrift/wiki/Graphical-User-Interface'
        print('Opening help website:\n' + help_url)
        import webbrowser
        webbrowser.open(help_url)

    def check_seeding(self):
        print('#'*50)
        print('Hang on, plot is comming in a few seconds...')
        mapres = self.mapresvar.get()[0]
        if mapres == 'f':
            print('...actually more like 30 seconds for full resolution coastline....')
            if self.has_diana is True:
                print('Du far ta deg ein liten trall mens du ventar.')
        print('#'*50)
        month = np.int(self.months.index(self.monthvar.get()) + 1)
        start_time = datetime(np.int(self.yearvar.get()), month,
                              np.int(self.datevar.get()),
                              np.int(self.hourvar.get()),
                              np.int(self.minutevar.get()))
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

        so = Leeway(loglevel=50)
        so.seed_elements(lon=lon, lat=lat, number=5000,
                         radius=radius, time=start_time)
        so.set_config('general:basemap_resolution', mapres)         
        so.plot(buffer=.5)
        del so

    def run_opendrift(self):
        sys.stdout.write('running OpenDrift')
        try:
            self.budgetbutton.destroy()
        except Exception as e:
            print(e)
            pass
        month = np.int(self.months.index(self.monthvar.get()) + 1)
        start_time = datetime(np.int(self.yearvar.get()), month,
                              np.int(self.datevar.get()),
                              np.int(self.hourvar.get()),
                              np.int(self.minutevar.get()))
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
                    print('Leeway object category: ' + lc)
                    break
            o.seed_elements(lon=lon, lat=lat, number=5000,
                            radius=radius, time=start_time,
                            objectType=ln + 1)
        if self.model.get() == 'OpenOil':
            o = OpenOil3D(weathering_model='noaa', loglevel=0)
            o.seed_elements(lon=lon, lat=lat, number=5000, radius=radius,
                            time=start_time, cone=cone,
                            oiltype=self.oljetype.get())

        readers = [  # Note that order (priority) is important!
            '/lustre/storeA/project/copernicus/sea/romsnorkyst/zdepths1h/*fc*.nc',
            'http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be',
            '/lustre/storeA/project/copernicus/sea/romsnordic/zdepths1h/roms_nordic4_ZDEPTHS_hr.fc.*.nc',
            'http://thredds.met.no/thredds/dodsC/sea/nordic4km/zdepths1h/aggregate_be',
            '/lustre/storeA/project/metproduction/products/meps/symlinks/thredds/meps_det_pp_2_5km_latest.nc',
            'http://thredds.met.no/thredds/dodsC/meps25files/meps_det_pp_2_5km_latest.nc',
            '/lustre/storeA/project/metproduction/products/arome2_5_arctic/thredds/arome_arctic_pp_2_5km_latest.nc',
            'http://thredds.met.no/thredds/dodsC/aromearcticlatest/arome_arctic_pp_2_5km_latest.nc',
            '/lustre/storeA/project/copernicus/sea/mywavewam4/*fc*.nc',
            'http://thredds.met.no/thredds/dodsC/sea/mywavewam4/mywavewam4_be',
            'http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_91.2/uv3z',
            'http://oos.soest.hawaii.edu/thredds/dodsC/hioos/model/atm/ncep_global/NCEP_Global_Atmospheric_Model_best.ncd']
        o.add_readers_from_list(readers)
        o.set_config('general:basemap_resolution', 'h')

        time_step = 1800  # Half hour
        duration = int(self.durationhours.get())*3600/time_step
        if self.directionvar.get() == 'backwards':
            time_step = -time_step
        if self.has_diana is True:
            extra_args = {'outfile': self.outputdir + '/opendrift_' +                          self.model.get() + o.start_time.strftime(
                                '_%Y%m%d_%H%M.nc')}
        else:
            extra_args = {}

        mapres = self.mapresvar.get()[0]
        o.set_config('general:basemap_resolution', mapres)         

        o.run(steps=duration, time_step=time_step,
              time_step_output=time_step, **extra_args)
        print(o)

        if self.has_diana is True:
            diana_filename = self.dianadir + '/opendrift_' + \
                self.model.get() + o.start_time.strftime(
                                '_%Y%m%d_%H%M.nc')
            o.write_netcdf_density_map(diana_filename)
            tk.Button(self.seed, text='Show in Diana',
                      command=lambda: os.system('diana &')
                      ).grid(row=8, column=2, sticky=tk.W, pady=4)
        tk.Button(self.seed, text='Animation',
                  command=o.animation).grid(row=8, column=3,
                                            sticky=tk.W, pady=4)
        if self.model.get() == 'OpenOil':
            self.budgetbutton = tk.Button(self.seed,
                text='Oil Budget', command=o.plot_oil_budget)
            self.budgetbutton.grid(row=8, column=4, sticky=tk.W, pady=4)


if __name__ == '__main__':
    OpenDriftGUI().mainloop()
