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
from importlib import resources
import opendrift
from opendrift.models.oceandrift import OceanDrift
from opendrift.models.openoil import OpenOil
from opendrift.models.leeway import Leeway
from opendrift.models.shipdrift import ShipDrift
from opendrift.models.openberg_old import OpenBergOld
from opendrift.models.plastdrift import PlastDrift
from opendrift.models.radionuclides import RadionuclideDrift

# Class to redirect output to text box
class TextRedirector:

    def __init__(self, widget, tag='stdout'):
        self.defstdout = sys.stdout
        self.widget = widget
        self.tag = tag

    def write(self, str):
        self.widget.configure(state='normal')
        self.widget.insert('end', str, (self.tag,))
        self.widget.update_idletasks()
        self.widget.see(tk.END)

    def flush(self):
        self.defstdout.flush()

# Class for help-text
# https://stackoverflow.com/questions/20399243/display-message-when-hovering-over-something-with-mouse-cursor-in-python
class ToolTip(object):

    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    def showtip(self, text):
        "Display text in tooltip window"
        self.text = text
        if self.tipwindow or not self.text:
            return
        x, y, cx, cy = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 57
        y = y + cy + self.widget.winfo_rooty() +27
        self.tipwindow = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        label = tk.Label(tw, text=self.text, justify=tk.LEFT,
                      background="#ffff00", relief=tk.SOLID, borderwidth=1,
                      font=("tahoma", "10", "normal"))
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

def CreateToolTip(widget, text):
    toolTip = ToolTip(widget)
    def enter(event):
        toolTip.showtip(text)
    def leave(event):
        toolTip.hidetip()
    widget.bind('<Enter>', enter)
    widget.bind('<Leave>', leave)


class OpenDriftGUI(tk.Tk):

    # Supported models as dictionary {model_name:model_class}
    opendrift_models = {m.__name__:m for m in
        [Leeway, OpenOil, ShipDrift, OpenBergOld, OceanDrift, PlastDrift, RadionuclideDrift]}

    extra_args = {'OpenOil': {'location': 'NORWAY'}}

    # Overriding some default config settings, suitable for GUI
    # TODO: should be set as default-default
    GUI_config = {
            'general:time_step_minutes': {'default': 15, 'min': 1},
            'general:time_step_output_minutes': {'default': 30, 'min': 5},
            'seed:number': {'default': 5000, 'max': 100000},
            'seed:m3_per_hour': {'default': 100}
            }

    def __init__(self):

        tk.Tk.__init__(self)

        self.title('OpenDrift ' + opendrift.__version__ + ' GTI Turbo Ultra')

        ##################
        # Layout frames
        ##################
        self.n = ttk.Notebook(self.master)
        self.n.grid()
        self.seed = ttk.Frame(self.n)
        self.confignotebook = ttk.Notebook(self.n)
        self.config = ttk.Frame(self.confignotebook)
        self.forcing = ttk.Frame(self.n)
        self.n.add(self.seed, text='Seeding')
        self.n.add(self.confignotebook, text='Config')
        self.n.add(self.forcing, text='Forcing')
        self.confignotebook.add(self.config, text='SubConfig')

        # Top
        self.top = tk.Frame(self.seed,
                            relief=tk.FLAT, pady=25, padx=25)
        self.top.grid(row=0, column=1, rowspan=1)
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

        #######################################################
        tk.Label(self.top, text='Simulation type').grid(row=0, column=0)
        self.model = tk.StringVar()
        self.model.set(list(self.opendrift_models)[0])
        self.modeldrop = tk.OptionMenu(self.top, self.model,
            *(list(self.opendrift_models)), command=self.set_model)
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
                            width=10, justify=tk.RIGHT)
        self.lon = tk.Entry(self.start, textvariable=self.lonvar,
                            width=10, justify=tk.RIGHT)
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
        conv=tk.Label(self.start, text='Convert from deg/min/sec', fg='blue')
        conv.grid(row=11, column=0, columnspan=2)
        conv.bind("<Button-1>", lambda e: self.convert_lonlat())

        ##########
        # Time
        ##########
        now = datetime.utcnow()
        tk.Label(self.start, text='Day').grid(row=20, column=0)
        tk.Label(self.start, text='Month').grid(row=20, column=1)
        tk.Label(self.start, text='Year').grid(row=20, column=2)
        tk.Label(self.start, text='Hour').grid(row=20, column=3)
        tk.Label(self.start, text='Minutes').grid(row=20, column=4)
        tk.Label(self.start, text='Timezone').grid(row=0, column=4)
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
        self.years = range(2015, now.year+2)
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

        self.timezonevar = tk.StringVar()
        self.timezone = ['UTC', 'CET']
        self.timezonevar.set('UTC')
        self.timezone = tk.OptionMenu(self.start, self.timezonevar,
                                    *self.timezone)
        self.timezone.grid(row=10, column=4)

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
        self.elat = tk.Entry(self.end, width=10, justify=tk.RIGHT)
        self.elon = tk.Entry(self.end, width=10, justify=tk.RIGHT)
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
        tk.Label(self.end, text='Minutes', bg='gray').grid(row=20, column=4)
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
        self.eyears = range(2015, now.year+2)
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
        if os.getenv('OPENDRIFT_GUI_OUTPUT', 'gui') == 'gui':
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
        self.set_model(list(self.opendrift_models)[0])

        with resources.open_text('opendrift.scripts', 'data_sources.txt') as fd:
            forcingfiles = fd.readlines()

        print(forcingfiles)
        for i, ff in enumerate(forcingfiles):
            tk.Label(self.forcing, text=ff.strip(), wraplength=650, font=('Courier', 8)).grid(
                     row=i, column=0, sticky=tk.W)

        ##########################
        try:
            img = ImageTk.PhotoImage(Image.open(
                self.o.test_data_folder() +
                                     '../../docs/opendrift_logo.png'))
            panel = tk.Label(self.seed, image=img)
            panel.image = img
            panel.grid(row=0, column=0)
        except Exception as e:
            print(e)
            pass # Could not display logo

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

        elif command == 'showanimationspecie':
            self.o.animation(
                color='specie',vmin=0,vmax=self.o.nspecies-1,
                colorbar=False,
                legend=[self.o.specie_num2name(i) for i in range(self.o.nspecies)]
                )
        elif command == 'saveanimationspecie':
            filename = filename + '.mp4'
            self.o.animation(filename=filename,
                color='specie',vmin=0,vmax=self.o.nspecies-1,
                                colorbar=False,
                legend=[self.o.specie_num2name(i) for i in range(self.o.nspecies)]
                )
        elif command == 'saveconcfile':
            self.o.guipp_saveconcfile(filename=homefolder+'/conc_radio.nc')
        elif command == 'plotconc':
            zlayer = [-1]
            time   = None
            specie = ['LMM']
            self.o.guipp_plotandsaveconc(filename=homefolder+'/conc_radio.nc',
                                         outfilename=homefolder+'/radio_plots/RadioConc',
                                         zlayers=zlayer, time=time, specie=specie )
        elif command == 'showanimationprofile':
            self.o.guipp_showanimationprofile()
        elif command == 'copy_netcdf':
            import shutil
            from tkinter import filedialog
            import pathlib
            folder_selected = filedialog.askdirectory(initialdir=pathlib.Path.home())
            try:
                shutil.copy(self.o.outfile_name, folder_selected)
                print('Copied netCDF file to:\n' \
                      f'{folder_selected}/{os.path.basename(self.o.outfile_name)}')
            except Exception as e:
                print('Could not copy file:')
                print(e)

    def validate_config(self, value_if_allowed, prior_value, key):
        """From config menu selection."""
        print('Input: %s -> %s', (key, value_if_allowed))
        if value_if_allowed == 'None':
            print('Setting None value')
            return True
        if value_if_allowed in ' -':
            print('Allowing temporally empty or minus sign')
            return True
        sc = self.o._config[key]
        if sc['type'] in ['int', 'float']:
            try:
                value_if_allowed = float(value_if_allowed)
            except:
                print('Nonumber')
                return False
        try:
            print('Setting: %s -> %s', (key, value_if_allowed))
            self.o.set_config(key, value_if_allowed)
            return True
        except:
            return False

    def set_model(self, model, rebuild_gui=True):

        # Creating simulation object (self.o) of chosen model class
        print('Setting model: ' + model)
        if model in self.extra_args:
            extra_args = self.extra_args[model]
        else:
            extra_args = {}
        self.o = self.opendrift_models[model](**extra_args)
        self.modelname = model  # So that new instance may be initiated at repeated run

        # Setting GUI-specific default config values
        for k,v in self.GUI_config.items():
            try:
                if 'default' in v:
                    self.o._set_config_default(k, v['default'])
                if 'min' in v:
                    self.o._config[k]['min'] = v['min']
                if 'max' in v:
                    self.o._config[k]['max'] = v['max']
            except:
                pass

        if rebuild_gui is False:
            return

        # Remove current GUI components and rebuild with new
        for con in self.confignotebook.winfo_children():
            con.destroy()
        self.subconfig = {}
        confnames = list(set([cn.split(':')[0] for cn in self.o._config]))
        confnames.extend(['environment:constant', 'environment:fallback'])
        confnames.remove('environment')
        for sub in confnames:
            self.subconfig[sub] = tk.Frame(self.confignotebook, pady=25)
            self.confignotebook.add(self.subconfig[sub], text=sub)

        sc = self.o.get_configspec(level=[2, 3])
        self.config_input = {}
        self.config_input_var = {}
        for i, key in enumerate(list(sc)):
            if key.startswith('environment:constant'):
                tab = self.subconfig['environment:constant']
                keystr = key.split(':')[-1]
            elif key.startswith('environment:fallback'):
                tab = self.subconfig['environment:fallback']
                keystr = key.split(':')[-1]
            else:
                tab = self.subconfig[key.split(':')[0]]
                keystr = ''.join(key.split(':')[1:])
            if keystr == '':
                keystr = key
            lab = tk.Label(tab, text=keystr)
            lab.grid(row=i, column=1, rowspan=1)
            if sc[key]['type'] in ['float', 'int']:
                self.config_input_var[i] = tk.StringVar()
                vcmd = (tab.register(self.validate_config),
                    '%P', '%s', key)
                self.config_input[i] = tk.Entry(
                    tab, textvariable=self.config_input_var[i],
                    validate='key', validatecommand=vcmd,
                    width=12, justify=tk.RIGHT)
                self.config_input[i].insert(0, str(sc[key]['default']))
                self.config_input[i].grid(row=i, column=2, rowspan=1)
                tk.Label(tab, text='[%s]  min: %s, max: %s' % (
                    sc[key]['units'], sc[key]['min'], sc[key]['max'])
                        ).grid(row=i, column=3, rowspan=1)
            if sc[key]['type'] == 'str':
                self.config_input_var[i] = tk.StringVar()
                vcmd = (tab.register(self.validate_config),
                    '%P', '%s', key)
                max_length = sc[key].get('max_length') or 12
                max_length = np.minimum(max_length, 64)
                self.config_input[i] = tk.Entry(
                    tab, textvariable=self.config_input_var[i],
                    validate='key', validatecommand=vcmd,
                    width=max_length, justify=tk.RIGHT)
                self.config_input[i].insert(0, str(sc[key]['default']))
                self.config_input[i].grid(row=i, column=2, columnspan=3, rowspan=1)
                #tk.Label(tab, text='').grid(row=i, column=3, rowspan=1)
            elif sc[key]['type'] == 'bool':
                if self.o.get_config(key) is True:
                    value = 1
                else:
                    value = 0
                self.config_input_var[i] = tk.IntVar(value=value)
                vcb = (tab.register(self.set_config_checkbox),
                       key, i)
                self.config_input[i] = tk.Checkbutton(
                    tab, variable=self.config_input_var[i],
                    command=vcb, text='')
                self.config_input[i].grid(row=i, column=2, rowspan=1)
            elif sc[key]['type'] == 'enum':
                self.config_input_var[i] = tk.StringVar(value=self.o.get_config(key))
                width = len(max(sc[key]['enum'], key=len))
                self.config_input[i] = ttk.Combobox(
                    tab, width=width,
                    textvariable=self.config_input_var[i],
                    values=sc[key]['enum'])
                self.config_input[i].bind("<<ComboboxSelected>>",
                        lambda event, keyx=key, ix=i:
                            self.set_config_enum(event, keyx, ix))
                self.config_input[i].grid(row=i, column=2, rowspan=1)

            CreateToolTip(lab, sc[key]['description'])

        try:
            self.results.destroy()
        except:
            pass

        # Only ESSENTIAL config items are shown on front page with seeding
        sc = self.o.get_configspec(level=opendrift.config.CONFIG_LEVEL_ESSENTIAL)
        self.seed_input = {}
        self.seed_input_var = {}
        self.seed_input_label = {}
        self.seed_frame = tk.Frame(self.seed, bd=2,
                                   relief=tk.FLAT, padx=5, pady=0)
        self.seed_frame.grid(row=60, columnspan=8, sticky='nsew')
        # FIND
        for num, i in enumerate(sc):
            varlabel = i.split(':')[-1]
            if i in self.o.ElementType.variables.keys():
                if 'units' in self.o.ElementType.variables[i].keys():
                    units = self.o.ElementType.variables[i]['units']
                    if units == '1':
                        units = 'fraction'
                    varlabel = '%s [%s]' % (varlabel, units)

            self.seed_input_label[i] = tk.Label(self.seed_frame,
                                                text=varlabel + '\t')
            self.seed_input_label[i].grid(row=num, column=0)
            CreateToolTip(self.seed_input_label[i], text=sc[i]['description'])
            actual_val = self.o.get_config(i)
            if sc[i]['type'] == 'enum':
                self.seed_input_var[i] = tk.StringVar()
                self.seed_input[i] = ttk.Combobox(
                    self.seed_frame, width=50,
                    textvariable=self.seed_input_var[i],
                    values=sc[i]['enum'])
                self.seed_input_var[i].set(actual_val)
            elif sc[i]['type'] == 'bool':
                self.seed_input_var[i] = tk.IntVar(value=sc[i]['value'])
                self.seed_input[i] = tk.Checkbutton(
                    self.seed_frame, variable=self.seed_input_var[i],
                    text=sc[i]['description'])
            else:
                self.seed_input_var[i] = tk.StringVar()
                max_length = sc[i].get('max_length') or 12
                max_length = np.minimum(max_length, 64)
                self.seed_input[i] = tk.Entry(
                    self.seed_frame, textvariable=self.seed_input_var[i],
                    width=max_length, justify=tk.RIGHT)
                self.seed_input[i].insert(0, actual_val)
            self.seed_input[i].grid(row=num, column=1)

    def set_config_checkbox(self, key, i):
        i = int(i)
        newval = self.config_input_var[i].get()
        if newval == 0:
            print('Setting %s to False' % key)
            self.o.set_config(key, False)
        elif newval == 1:
            print('Setting %s to True' % key)
            self.o.set_config(key, True)

    def set_config_enum(self, event, key, i):
        newval = self.config_input_var[i].get()
        print('Setting ' + key + newval)
        self.o.set_config(key, newval)

    def show_help(self):
        help_url = 'https://opendrift.github.io/gui.html'
        print('Opening help website:\n' + help_url)
        import webbrowser
        webbrowser.open(help_url)

    def convert_lonlat(self):
        convert_url = 'https://www.rapidtables.com/convert/number/degrees-minutes-seconds-to-degrees.html'
        print('Opening conversion website:\n' + convert_url)
        import webbrowser
        webbrowser.open(convert_url)

    def check_seeding(self):
        print('#'*50)
        print('Hang on, plot is comming in a few seconds...')
        print('#'*50)
        month = int(self.months.index(self.monthvar.get()) + 1)
        start_time = datetime(int(self.yearvar.get()), month,
                              int(self.datevar.get()),
                              int(self.hourvar.get()),
                              int(self.minutevar.get()))
        emonth = int(self.months.index(self.emonthvar.get()) + 1)
        end_time = datetime(int(self.eyearvar.get()), emonth,
                            int(self.edatevar.get()),
                            int(self.ehourvar.get()),
                            int(self.eminutevar.get()))
        sys.stdout.flush()
        lon = float(self.lon.get())
        lat = float(self.lat.get())
        radius = float(self.radius.get())
        elon = float(self.elon.get())
        elat = float(self.elat.get())
        eradius = float(self.eradius.get())
        if lon != elon or lat != elat or start_time != end_time:
            lon = [lon, elon]
            lat = [lat, elat]
            radius = [radius, eradius]
            start_time = [start_time, end_time]
            cone = True
        else:
            cone = False

        so = Leeway(loglevel=50)
        for k,v in self.GUI_config.items():
            try:
                so.set_config(k, v)
            except:
                pass
        so.seed_cone(lon=lon, lat=lat, radius=radius, time=start_time)
        so.plot(buffer=.5, fast=True)
        del so

    def run_opendrift(self):
        sys.stdout.write('running OpenDrift')

        try:
            self.budgetbutton.destroy()
        except Exception as e:
            #print(e)
            pass
        month = int(self.months.index(self.monthvar.get()) + 1)
        start_time = datetime(int(self.yearvar.get()), month,
                              int(self.datevar.get()),
                              int(self.hourvar.get()),
                              int(self.minutevar.get()))
        emonth = int(self.months.index(self.emonthvar.get()) + 1)
        end_time = datetime(int(self.eyearvar.get()), emonth,
                            int(self.edatevar.get()),
                            int(self.ehourvar.get()),
                            int(self.eminutevar.get()))

        timezone = self.timezonevar.get()
        if timezone != 'UTC':
            import pytz
            local = pytz.timezone(timezone)
            local_start = local.localize(start_time, is_dst=None)
            local_end = local.localize(end_time, is_dst=None)
            start_time = local_start.astimezone(pytz.utc).replace(tzinfo=None)
            end_time = local_end.astimezone(pytz.utc).replace(tzinfo=None)

        # Creating fresh instance of the current model, but keeping config
        adjusted_config = self.o._config
        self.set_model(self.modelname, rebuild_gui=False)
        self.o._config = adjusted_config

        # Add secondary log-handler to file
        if self.has_diana is True and True:
            logfile = self.outputdir + '/opendrift_' + self.modelname + start_time.strftime('_%Y%m%d_%H%M.log')
            import logging
            if hasattr(self, 'handler'):
                logging.getLogger('opendrift').removeHandler(self.handler)
            self.handler = logging.FileHandler(logfile, mode='w')
            self.handler.setFormatter(logging.getLogger('opendrift').handlers[0].formatter)
            logging.getLogger('opendrift').addHandler(self.handler)
            if len(logging.getLogger('opendrift').handlers) > 2:
                raise ValueError('Too many log handlers')

        sys.stdout.flush()
        lon = float(self.lon.get())
        lat = float(self.lat.get())
        radius = float(self.radius.get())
        elon = float(self.elon.get())
        elat = float(self.elat.get())
        eradius = float(self.eradius.get())
        if lon != elon or lat != elat or start_time != end_time:
            lon = [lon, elon]
            lat = [lat, elat]
            radius = [radius, eradius]
            start_time = [start_time, end_time]
            cone = True
        else:
            cone = False

        for se in self.seed_input:
            val = self.seed_input_var[se].get()
            if self.o._config[se]['type'] in ['float', 'int']:
                val = float(val)
            elif self.o._config[se]['type'] == 'bool':
                if val == 1:
                    val = True
                elif val == 0:
                    val = False
                else:
                    nothing
            self.o.set_config(se, val)

        with resources.path('opendrift.scripts', 'data_sources.txt') as f:
            self.o.add_readers_from_file(f)

        self.o.seed_cone(lon=lon, lat=lat, radius=radius,
                         time=start_time)#, #cone=cone,
                         #**extra_seed_args)

        time_step = self.o.get_config('general:time_step_minutes')*60
        time_step_output = self.o.get_config('general:time_step_output_minutes')*60
        duration = int(self.durationhours.get())*3600/time_step
        extra_args = {'time_step': time_step, 'time_step_output': time_step_output}
        if self.directionvar.get() == 'backwards':
            extra_args['time_step'] = -extra_args['time_step']
            extra_args['time_step_output'] = -extra_args['time_step_output']
        if self.has_diana is True:
            extra_args['outfile'] = self.outputdir + '/opendrift_' + \
                self.model.get() + self.o.start_time.strftime('_%Y%m%d_%H%M.nc')

        self.simulationname = 'opendrift_' + self.model.get() + \
            self.o.start_time.strftime('_%Y%m%d_%H%M')

        # Starting simulation run
        self.o.run(steps=duration, **extra_args)
        print(self.o)

        try:
            os.chmod(extra_args['outfile'], 0o666)
        except:
            pass

        # Model-specific post processing
        self.o.gui_postproc()


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


        if self.model.get() =='RadionuclideDrift':
            tk.Button(self.results, text='Show animation specie',
                      command=lambda: self.handle_result(
                          'showanimationspecie')).grid(row=30, column=2)
            tk.Button(self.results, text='Save animation specie',
                      command=lambda: self.handle_result(
                          'saveanimationspecie')).grid(row=40, column=2)
            tk.Button(self.results, text='Animation profile',
                      command=lambda: self.handle_result(
                          'showanimationprofile')).grid(row=10, column=2)
#             tk.Button(self.results, text='Save conc file',
#                       command=lambda: self.handle_result(
#                           'saveconcfile')).grid(row=20, column=2)
            tk.Button(self.results, text='Plot conc',
                      command=lambda: self.handle_result(
                          'plotconc')).grid(row=20, column=2)



        if self.has_diana is True:
            diana_filename = self.dianadir + self.simulationname + '.nc'
            self.o.write_netcdf_density_map(diana_filename)
            tk.Button(self.results, text='Show in Diana',
                      command=lambda: os.system('diana &')
                      ).grid(row=80, column=1)

            try:
                os.chmod(diana_filename, 0o666)
            except:
                pass

            tk.Button(self.results, text='Copy netCDF file',
                      command=lambda: self.handle_result(
                          'copy_netcdf')).grid(row=81, column=1)

def main():
    OpenDriftGUI().mainloop()

if __name__ == '__main__':

    # Add any argument to redirect output to terminal instead of GUI window.
    # TODO: must be a better way to pass arguments to Tkinter?
    if len(sys.argv) > 1:
        os.environ['OPENDRIFT_GUI_OUTPUT'] = 'terminal'
    else:
        os.environ['OPENDRIFT_GUI_OUTPUT'] = 'gui'
    OpenDriftGUI().mainloop()
