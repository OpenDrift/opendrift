#!/usr/bin/env python

import matplotlib
from matplotlib import pyplot as plt
if __name__ == '__main__':
    matplotlib.use('TKAgg')

import sys
import os
from datetime import datetime, timedelta
import numpy as np

from PIL import ImageTk, Image
import tkinter as tk
from tkinter import ttk
import opendrift
from opendrift.models.openoil3D import OpenOil3D
from opendrift.models.leeway import Leeway
from opendrift.models.shipdrift import ShipDrift


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
        self.start_t.grid(row=20, column=0, rowspan=1)
        self.end_t = tk.Frame(self.seed, relief=tk.FLAT)
        self.end_t.grid(row=30, column=0, rowspan=1)
        self.start = tk.Frame(self.seed, bg='lightgray', bd=2,
                              relief=tk.SUNKEN, pady=5, padx=5)
        self.start.grid(row=20, column=1, rowspan=1)
        self.end = tk.Frame(self.seed, bg='gray', bd=2,
                            relief=tk.SUNKEN, padx=5, pady=5)
        self.end.grid(row=30, column=1)
        self.coastline = tk.Frame(self.seed, bd=2,
                                 relief=tk.FLAT, padx=5, pady=0)
        self.coastline.grid(row=40, column=1)
        self.duration = tk.Frame(self.seed, bd=2,
                                 relief=tk.FLAT, padx=5, pady=5)
        self.duration.grid(row=50, column=1)
        self.output = tk.Frame(self.seed, bd=2,
                               relief=tk.FLAT, padx=5, pady=0)
        self.output.grid(row=70, column=0, columnspan=7, sticky='nsew')

        self.results = tk.Frame(self.seed, bd=2,
                               relief=tk.FLAT, padx=5, pady=0)
        self.results.grid(row=60, column=7, columnspan=1, sticky='ew')

        ##########################
        self.title('OpenDrift')
        self.o = OpenOil3D(weathering_model='noaa', location='NORWAY')
        try:
            img = ImageTk.PhotoImage(Image.open(self.o.test_data_folder() +
                                     '../../docs/opendrift_logo.png'))
            panel = tk.Label(self.seed, image=img)
            panel.image = img
            panel.grid(row=0, column=0)
        except:
            pass # Could not display logo
        #######################################################
        tk.Label(self.top, text='Simulation type').grid(row=0, column=0)
        self.model = tk.StringVar()
        models = opendrift.get_model_names()
        self.model.set(models[0])
        self.modeldrop = tk.OptionMenu(self.top, self.model,
                                       *(models), command=self.set_model)
        self.modeldrop.grid(row=0, column=1)

        help_button = tk.Button(self.top, text='Help',
                                command=self.show_help)
        help_button.grid(row=0, column=2, padx=50)

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
        self.lon.grid(row=10, column=1)
        self.lon.insert(0, '4.5')
        self.lat.grid(row=10, column=0)
        self.lat.insert(0, '60.0')
        self.radius.grid(row=10, column=2)
        self.radius.insert(0, '1000')
        self.lonvar.trace('w', self.copy_position)
        self.latvar.trace('w', self.copy_position)
        self.radiusvar.trace('w', self.copy_position)
        ##########
        # Time
        ##########
        now = datetime.utcnow()
        tk.Label(self.start, text='Day').grid(row=20, column=0)
        tk.Label(self.start, text='Month').grid(row=20, column=1)
        tk.Label(self.start, text='Year').grid(row=20, column=2)
        tk.Label(self.start, text='Hour').grid(row=20, column=3)
        tk.Label(self.start, text='Minutes [UTC]').grid(row=20, column=4)
        self.datevar = tk.StringVar()
        self.dates = range(1, 32)
        self.datevar.set(now.day)
        self.date = tk.OptionMenu(self.start, self.datevar, *self.dates)
        self.date.grid(row=30, column=0)

        self.monthvar = tk.StringVar()
        self.months = ['January', 'February', 'March', 'April', 'May',
                       'June', 'July', 'August', 'September', 'October',
                       'November', 'December']
        self.monthvar.set(self.months[now.month-1])
        self.month = tk.OptionMenu(self.start, self.monthvar,
                                   *self.months)
        self.month.grid(row=30, column=1)

        self.yearvar = tk.StringVar()
        self.years = range(2015, now.year+1)
        self.yearvar.set(now.year)
        self.year = tk.OptionMenu(self.start, self.yearvar, *self.years)
        self.year.grid(row=30, column=2)

        self.hourvar = tk.StringVar()
        self.hours = range(0, 24)
        self.hourvar.set(now.hour)
        self.hour = tk.OptionMenu(self.start, self.hourvar, *self.hours)
        self.hour.grid(row=30, column=3)

        self.minutevar = tk.StringVar()
        self.minutes = range(0, 60, 5)
        self.minutevar.set(now.minute)
        self.minute = tk.OptionMenu(self.start, self.minutevar,
                                    *self.minutes)
        self.minute.grid(row=30, column=4)

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
        self.elon.grid(row=10, column=1)
        self.elon.insert(0, '4.5')
        self.elat.grid(row=10, column=0)
        self.elat.insert(0, '60.0')
        self.eradius.grid(row=10, column=2)
        self.eradius.insert(0, '1000')
        ##########
        # Time
        ##########
        now = datetime.utcnow()
        tk.Label(self.end, text='Day', bg='gray').grid(row=20, column=0)
        tk.Label(self.end, text='Month', bg='gray').grid(row=20, column=1)
        tk.Label(self.end, text='Year', bg='gray').grid(row=20, column=2)
        tk.Label(self.end, text='Hour', bg='gray').grid(row=20, column=3)
        tk.Label(self.end, text='Minutes [UTC]', bg='gray').grid(row=20, column=4)
        self.edatevar = tk.StringVar()
        self.edates = range(1, 32)
        self.edatevar.set(now.day)
        self.edate = tk.OptionMenu(self.end, self.edatevar, *self.edates)
        self.edate.grid(row=30, column=0)

        self.emonthvar = tk.StringVar()
        self.emonthvar.set(self.months[now.month-1])
        self.emonth = tk.OptionMenu(self.end, self.emonthvar,
                                    *self.months)
        self.emonth.grid(row=30, column=1)

        self.eyearvar = tk.StringVar()
        self.eyears = range(2015, now.year+1)
        self.eyearvar.set(now.year)
        self.eyear = tk.OptionMenu(self.end, self.eyearvar, *self.eyears)
        self.eyear.grid(row=30, column=2)

        self.ehourvar = tk.StringVar()
        self.ehours = range(0, 24)
        self.ehourvar.set(now.hour)
        self.ehour = tk.OptionMenu(self.end, self.ehourvar, *self.ehours)
        self.ehour.grid(row=30, column=3)

        self.eminutevar = tk.StringVar()
        self.eminutes = range(0, 60, 5)
        self.eminutevar.set(now.minute)
        self.eminute = tk.OptionMenu(self.end, self.eminutevar,
                                     *self.eminutes)
        self.eminute.grid(row=30, column=4)
        self.eyear.config(bg='gray')
        self.emonth.config(bg='gray')
        self.edate.config(bg='gray')
        self.ehour.config(bg='gray')
        self.eminute.config(bg='gray')

        # Check seeding
        check_seed = tk.Button(self.end_t, text='Check seeding',
                                command=self.check_seeding)
        check_seed.grid(row=10, column=0, padx=0)

        #######################
        # Simulation duration
        #######################
        #tk.Label(self.coastline, text='Coastline resolution ').grid(
        #         row=40, column=1)
        #self.mapresvar = tk.StringVar()
        #self.mapres = tk.OptionMenu(self.coastline, self.mapresvar,
        #                            *['full', 'high'])
        #self.mapres.grid(row=40, column=2)
        #self.mapresvar.set('high')

        tk.Label(self.duration, text='Run simulation ').grid(row=50, column=0)
        self.durationhours = tk.Entry(self.duration, width=3,
                                      justify=tk.RIGHT)
        self.durationhours.grid(row=50, column=1)
        self.durationhours.insert(0, 12)
        tk.Label(self.duration, text=' hours ').grid(row=50, column=2)

        self.directionvar = tk.StringVar()
        self.directionvar.set('forwards')
        self.direction = tk.OptionMenu(self.duration, self.directionvar,
                                       'forwards', 'backwards')
        self.direction.grid(row=50, column=3)
        tk.Label(self.duration, text=' in time ').grid(row=50, column=4)

        ##############
        # Output box
        ##############
        self.text = tk.Text(self.output, wrap="word", height=18)
        self.text.grid(row=60, columnspan=6, sticky='nsw')
        self.text.tag_configure("stderr", foreground="#b22222")
        sys.stdout = TextRedirector(self.text, "stdout")
        sys.stderr = TextRedirector(self.text, "stderr")
        s = tk.Scrollbar(self)
        s.grid(row=60, column=6, sticky='ns')
        s.config(command=self.text.yview)
        self.text.config(yscrollcommand=s.set)

        # Diana
        self.dianadir = '/vol/vvfelles/opendrift/output/'
        if os.path.exists(self.dianadir):
            self.has_diana = True
            print('Diana is available!')
            self.outputdir = '/vol/vvfelles/opendrift/output_native/'
            startbutton = 'PEIS PAO'
        else:
            self.has_diana = False
            startbutton = 'START'

        ##############
        # Initialise
        ##############
        #o = OpenOil3D()
        self.set_model(self.available_models[0])

        ##########
        # RUN
        ##########
        tk.Button(self.seed, text=startbutton, bg='green',
                  command=self.run_opendrift).grid(row=80, column=1,
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

    def handle_result(self, command):

        from os.path import expanduser
        homefolder = expanduser("~")
        filename = homefolder + '/' + self.simulationname

        if command[0:4] == 'save':
            plt.switch_backend('agg')
        elif command[0:4] == 'show':
            plt.switch_backend('TkAgg')

        if command == 'saveanimation':
            filename = filename + '.mp4'
            self.o.animation(filename=filename)
            print('='*30 + '\nAnimation saved to file:\n'
                  + filename + '\n' + '='*30)
        elif command == 'showanimation':
            self.o.animation()
        elif command == 'saveplot':
            filename = filename + '.png'
            self.o.plot(filename=filename)
            print('='*30 + '\nPlot saved to file:\n'
                  + filename + '\n' + '='*30)
        elif command == 'showplot':
            self.o.plot()
        elif command == 'showoilbudget':
            self.o.plot_oil_budget()
        elif command == 'saveoilbudget':
            filename = filename + '_oilbudget.png'
            self.o.plot_oil_budget(filename=filename)
            print('='*30 + '\nPlot saved to file: '
                  + filename + '\n' + '='*30)

    def set_model(self, model):
        if model == 'OpenOil':
            self.o = OpenOil3D(weathering_model='noaa', location='NORWAY')
        elif model == 'Leeway':
            self.o = Leeway()
        elif model == 'ShipDrift':
            self.o = ShipDrift()

        for con in self.config.winfo_children():
            con.destroy()
        self.con = tk.Label(self.config, text="\n\nConfiguration\n\n")
        self.con.grid(row=0, column=1, rowspan=1)
        for i, cs in enumerate(self.o._config_hashstrings()):
            tk.Label(self.config, text=cs).grid(row=i, column=1, rowspan=1)
        try:
            self.results.destroy()
        except:
            pass
        try:
            # Removeing depth input boxes
            self.depthlabel.destroy()
            self.depth.destroy()
            self.seafloor.destroy()
            self.amount.destroy()
            self.amountlabel.destroy()
        except:
            pass

        print('Setting model: ' + model)
        print(self.o.list_configspec())
        sc = self.o.get_seed_config()
        print(sc)
        self.seed_input = {}
        self.seed_input_var = {}
        self.seed_input_label = {}
        self.seed_frame = tk.Frame(self.seed, bd=2,
                                   relief=tk.FLAT, padx=5, pady=0)
        self.seed_frame.grid(row=60, columnspan=8, sticky='nsew')
        # FIND
        for num, i in enumerate(sc):
            if i in ['ocean_only', 'jibeProbability']:
                continue  # workaround, should be avoided in future
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

        if model == 'OpenOil':  # User may be able to chose release depth
            self.depthlabel = tk.Label(self.duration, text='Spill depth [m]')
            self.depthlabel.grid(row=60, column=0)
            self.depthvar = tk.StringVar()
            self.depth = tk.Entry(self.duration, textvariable=self.depthvar,
                                  width=6, justify=tk.RIGHT)
            self.depth.grid(row=60, column=2)
            self.depth.insert(0, '0')
            self.seafloorvar = tk.IntVar()
            self.seafloor = tk.Checkbutton(self.duration, variable=self.seafloorvar,
                                           text='seafloor',
                                           command=self.seafloorbutton)
            self.seafloor.grid(row=60, column=3)

            self.amountvar = tk.StringVar()
            self.amount = tk.Entry(self.duration, textvariable=self.amountvar,
                                  width=6, justify=tk.RIGHT)
            self.amount.grid(row=60, column=4)
            self.amount.insert(0, '100')
            self.amountlabel = tk.Label(self.duration, text='m3 (per hour)')
            self.amountlabel.grid(row=60, column=5)

    def seafloorbutton(self):
        if self.seafloorvar.get() == 1:
            self.depth.config(state='disabled')
        else:
            self.depth.config(state='normal')

    def show_help(self):
        help_url = 'https://github.com/OpenDrift/opendrift/wiki/Graphical-User-Interface'
        print('Opening help website:\n' + help_url)
        import webbrowser
        webbrowser.open(help_url)

    def check_seeding(self):
        print('#'*50)
        print('Hang on, plot is comming in a few seconds...')
        #mapres = self.mapresvar.get()[0]
        #if mapres == 'f':
        #    print('...actually more like 30 seconds for full resolution coastline....')
        #    if self.has_diana is True:
        #        print('Du far ta deg ein liten trall mens du ventar.')
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
        so.plot(buffer=.5, fast=True)
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
            self.o = Leeway(loglevel=0)
            #for ln, lc in enumerate(self.leewaycategories):
            #    if self.oljetype.get() == lc.strip().replace('>', ''):
            #        print('Leeway object category: ' + lc)
            #        break
            #extra_seed_args = {'objectType': ln + 1}
        elif self.model.get() == 'OpenOil':
            self.o = OpenOil3D(weathering_model='noaa', loglevel=0)
            #extra_seed_args = {'oiltype': self.oljetype.get()}
        elif self.model.get() == 'ShipDrift':
            self.o = ShipDrift(loglevel=0)

        extra_seed_args = {}
        for se in self.seed_input:
            if se == 'ocean_only':
                continue  # To be fixed/removed
            val = self.seed_input[se].get()
            try:
                extra_seed_args[se] = np.float(val)
            except:
                extra_seed_args[se] = val
            self.o.set_config('seed:' + se, val)
        if 'object_type' in extra_seed_args:
            extra_seed_args['objectType'] = extra_seed_args['object_type']
            del extra_seed_args['object_type']

        self.o.add_readers_from_file(self.o.test_data_folder() +
            '../../opendrift/scripts/data_sources.txt')

        extra_seed_args = {}
        if self.model.get() == 'OpenOil':
            if self.seafloorvar.get() == 1:
                z = 'seafloor'
            else:
                z = -np.abs(np.float(self.depthvar.get()))  # ensure negative z
            extra_seed_args['z'] = z
            extra_seed_args['m3_per_hour'] = np.float(self.amountvar.get())
        self.o.seed_elements(lon=lon, lat=lat, number=5000, radius=radius,
                        time=start_time, cone=cone,
                        **extra_seed_args)

        #time_step = o.get_config('general:time_step_minutes')*60
        time_step = 900  # 15 minutes
        time_step_output = timedelta(minutes=30)
        duration = int(self.durationhours.get())*3600/time_step
        if self.directionvar.get() == 'backwards':
            time_step = -time_step
        if self.has_diana is True:
            extra_args = {'outfile': self.outputdir + '/opendrift_' +
                self.model.get() + self.o.start_time.strftime('_%Y%m%d_%H%M.nc')}
        else:
            extra_args = {}

        #mapres = self.mapresvar.get()[0]
        self.simulationname = 'opendrift_' + self.model.get() + \
            self.o.start_time.strftime('_%Y%m%d_%H%M')

        # Starting simulation run
        self.o.run(steps=duration, time_step=time_step,
              time_step_output=time_step_output, **extra_args)
        print(self.o)

        self.results.destroy()
        self.results = tk.Frame(self.seed, bd=2,
                               relief=tk.FLAT, padx=5, pady=0)
        self.results.grid(row=70, column=3, columnspan=1, sticky='ew')
        tk.Button(self.results, text='Show animation',
                  command=lambda: self.handle_result(
                    'showanimation')).grid(row=10, column=1)
        tk.Button(self.results, text='Save animation',
                  command=lambda: self.handle_result(
                    'saveanimation')).grid(row=20, column=1)
        tk.Button(self.results, text='Show plot',
                  command=lambda: self.handle_result(
                    'showplot')).grid(row=30, column=1)
        tk.Button(self.results, text='Save plot',
                  command=lambda: self.handle_result(
                    'saveplot')).grid(row=40, column=1)
        if self.model.get() == 'OpenOil':
            tk.Button(self.results, text='Save oil budget',
                      command=lambda: self.handle_result(
                        'saveoilbudget')).grid(row=50, column=1)
            tk.Button(self.results, text='Show oil budget',
                      command=lambda: self.handle_result(
                        'showoilbudget')).grid(row=60, column=1)

        if self.has_diana is True:
            diana_filename = self.dianadir + self.simulationname + '.nc'
            self.o.write_netcdf_density_map(diana_filename)
            tk.Button(self.results, text='Show in Diana',
                      command=lambda: os.system('diana &')
                      ).grid(row=80, column=1)

if __name__ == '__main__':
    OpenDriftGUI().mainloop()

