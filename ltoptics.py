import logging
import signal
import sys
import time
import datetime

import matplotlib.pyplot as plt
import numpy as np

from tkinter import *
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename

from lasertracker import *
from bam_LT_fxns import *
import bam_LT_fxns

import redis
import yaml
from functools import partial

# script 1
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

config_file = 'wfpt.yaml'

with open(config_file, 'r') as f:
    config = yaml.safe_load(f)

print(config)

lt_host = config['laser_tracker']['host']
lt_type = config['laser_tracker']['type']


if lt_type == 'leica_at4xx':
    import leica_at4xx as lt
else:
    print('Unknown laser tracker type', lt_type)
    exit(1)

# script 1
def signal_handler(signal, frame):
    logger.info('You pressed Ctrl+C! Exiting...')
    sys.exit(0)


# GUI commands
def fileSelect(*args):
    global running
    running = False

    global zemax_prescrip_fname
    zemax_prescrip_fname = askopenfilename(filetypes=[("txt files", ".txt")])

    Label(zemax_frame, text=(zemax_prescrip_fname), fg="Green").grid(columnspan=4, row=1, sticky=(W, E))


def ltlog(retro_id, zemax_surf = -1, mode = 'smr', coordsys='LT', xyz=[0,0,0], err=[0,0,0], zemax_file='', redis_connection=None):
    dt = datetime.datetime.now()
    dt_string = dt.strftime('%Y-%m-%dT%H:%M:%S')
    log_str = '%s\t%d\t%d\t%s\t%s\t%9.3f\t%9.3f\t%9.3f\t%6.3f\t%6.3f\t%6.3f\t%s\n' % ( dt_string, retro_id, zemax_surf, mode, coordsys, xyz[0], xyz[1], xyz[2],err[0],err[1],err[2],zemax_file)
    print('LOG: ', log_str)
    with open('lt.log', 'a') as f :
        f.write(log_str)
    if redis_connection is not None:
        redis_connection.set('lt.%s.%s.X' % (mode, coordsys), xyz[0])
        redis_connection.set('lt.%s.%s.Y' % (mode, coordsys), xyz[1])
        redis_connection.set('lt.%s.%s.Z' % (mode, coordsys), xyz[2])
  	
# measure button command, adds retro coordinates to cesapi_positions array
def measureRetro(laser_tracker, index, use_saved_position=False):
    global running
    running = False

    global cesapi_positions

    if use_saved_position:
        laser_tracker.goto_position(cesapi_positions[index])

    cesapi_positions[index] = laser_tracker.measure()
    global trackerpilot_positions
    trackerpilot_positions[index] = cesapi_positions[index]

    #autosave
    np.savetxt("tempbenchcoords.txt", cesapi_positions)

    ltlog(retro_id=index, mode='fixed', coordsys='LT', xyz=cesapi_positions[index], err=[laser_tracker.measurement.dStd1, laser_tracker.measurement.dStd2, laser_tracker.measurement.dStd3],redis_connection=myredis)

def measureAll(laser_tracker,max=None):
    if max==None:
        max = len(cesapi_positions)
    for i in range(max):
        if np.any(cesapi_positions[i] != [1, 1, 1]):
            print(i,cesapi_positions[i])
            measureRetro(laser_tracker, i,use_saved_position=True)


# SVD rotation to find unmeasured SMRs from measured
def SVD_Rotation(cesapi_positions, trackerpilot_positions, SMRs_measured_trackerpilot):

    print ("cesapi_positions: ")
    print(cesapi_positions)

    print ("trackerpilot_positions ")
    print(trackerpilot_positions)

    print('SMRs_measured_trackerpilot')
    print(SMRs_measured_trackerpilot)

    Q = np.ones((3,3))
    P = np.ones((3,3))
    SMRs_measured_trackerpilot_XYZ = np.ones((3,3))

    j = 0
    for i in range(len(cesapi_positions) - 1): 
        if np.any(cesapi_positions[i] != [1, 1, 1]):
            Q[j] = trackerpilot_positions[i]
            P[j] = SMRs_measured_trackerpilot[i]
            j = j+1

    for i in [0, 1, 2]:
        SMRs_measured_trackerpilot_XYZ[i] = SMRs_measured_trackerpilot[i]

    Q_vert = np.transpose(Q)
    P_vert = np.transpose(P)

    P_centroid = (P_vert[:, 0] + P_vert[:, 1] + P_vert[:, 2]) * (1 / 3)
    Q_centroid = (Q_vert[:, 0] + Q_vert[:, 1] + Q_vert[:, 2]) * (1 / 3)

    X = P_vert - P_centroid
    Y = Q_vert - Q_centroid

    S = X * np.transpose(Y)

    u, s, vh = np.linalg.svd(S)

    i = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, np.linalg.det(np.transpose(vh) * np.transpose(u))]])

    global R, translation
    global V
    V = SMRs_measured_trackerpilot_XYZ
    
    R = (np.transpose(vh) * i) * np.transpose(u)

    translation = Q_centroid - (R * P_centroid)
    new_q = (R * np.transpose(SMRs_measured_trackerpilot_XYZ)) + (translation)
    var = np.transpose(new_q)

    print (' ')
    print ('P, Q, new_q')
    print (P)
    print (Q)
    print (new_q)
    print (' ')
    print ('R t')
    print(R)
    print(translation)
    print(' ')
    print ('SMRs_measured_trackerpilot_XYZ')
    print(SMRs_measured_trackerpilot_XYZ)


    print(trackerpilot_positions)

    smr = trackerpilot_positions[5]

    print('Exit SVD_Rotation')

    return var, smr

# loads coordinates of smrs
def load():
    global running
    running = False

    global cesapi_positions
    global trackerpilot_positions

    file = askopenfilename(filetypes=[("txt files", ".txt")])

    cesapi_positions = np.loadtxt(file)

    trackerpilot_positions = np.loadtxt(file)

    print(cesapi_positions)
    print(trackerpilot_positions)


# saves coordinates of smrs
def save():
    global running
    running = False

    global cesapi_positions
    np.savetxt(asksaveasfilename(), cesapi_positions)


# resets cesapi_positions array
def reset():
    global running
    running = False

    global cesapi_positions
    cesapi_positions.fill(1.)
    print(cesapi_positions)

    global trackerpilot_positions
    trackerpilot_positions.fill(1.)
    print(trackerpilot_positions)

    np.savetxt("tempbenchcoords.txt", cesapi_positions)

#  adjusts mirror or SMR
def Adjust(mode="mirror", *args):

    print('In Adjust')
    
    global running

    global zemax_prescrip_fname
    global surface_number

    global nominal

    global smr_real_posn_LTXYZ
    global benchinfo
    global config

    print("In Adjust")

    # Save the current position
    # 
    place_to_measure = laser_tracker.measure()
    
    # Remeasure all the previously measured positions without human hands touching
    #
    print("Measuring previous positions")
    measureAll(laser_tracker,max=5)
    
    # Go back to the place we started
    #
    print("Back to where we started")
    laser_tracker.goto_position(place_to_measure)
    #command.GoPosition(True,place_to_measure.dVal1, place_to_measure.dVal2, place_to_measure.dVal3)

    data, smr = SVD_Rotation(cesapi_positions, trackerpilot_positions, permanent_smrs)

    if zemax_prescrip_fname == "No Zemax File Selected!":
        messagebox.showerror("ERROR", "No Zemax File Selected!")
    elif surface_number.get() == '':
        messagebox.showerror("ERROR", "No Surface Number Entered!")
    else:

        # Measure real (unreflected) SMR position
        smr_real_posn_LTXYZ = smr

        originsmrcoords_LTXYZ = np.array(data[0])[0]
        xaxissmrcoords_LTXYZ =  np.array(data[1])[0]
        yaxissmrcoords_LTXYZ =  np.array(data[2])[0]

        print(' ')
        print('Coordsys SMRs LTXYZ: ')
        print(originsmrcoords_LTXYZ)
        print(xaxissmrcoords_LTXYZ)
        print(yaxissmrcoords_LTXYZ)
        print(' ')
        
        print(config)
        
        Zemax_origin_in_benchcoordinates = np.array(config['Zemax_origin_in_benchcoordinates'])
        Zemax_unitvectors_in_benchcoordinates = np.array(config['Zemax_unitvectors_in_benchcoordinates'])

        benchinfo = bam_LT_fxns.measurement_set(originsmrcoords_LTXYZ, xaxissmrcoords_LTXYZ, yaxissmrcoords_LTXYZ,
                                                Zemax_origin_in_benchcoordinates,
                                                Zemax_unitvectors_in_benchcoordinates,
                                                zemax_prescrip_fname=zemax_prescrip_fname)

        LTcoords_bench_origin_smr = originsmrcoords_LTXYZ
        LTcoords_bench_xaxis_smr = xaxissmrcoords_LTXYZ
        LTcoords_bench_yaxis_smr = yaxissmrcoords_LTXYZ

        # Select temp.alignment flat
        benchinfo.select_surface(surface_number.get())

        benchinfo.measure_SMR_real(smr_real_posn_LTXYZ)
        print('--------------------------------------------------------------------------')
        benchinfo.get_virtualimg_posn()
        print('--------------------------------------------------------------------------')

        print('Mode = ', mode)
        if mode == "mirror":
            # Desired virtual image position:
            nominal = benchinfo.virtual_SMR_posn_LTcoord_Cart
        elif mode == "smr":
            nominal = benchinfo.transform_bench_to_LT(benchinfo.surface_coords_bench)

        print('Zemax origin in bench coordinates: ', Zemax_origin_in_benchcoordinates)
        print('Bench origin in Zemax coordinates: ', benchinfo.bench_origin_in_Zemaxcoordinates)
        running = True
        print('Running: ',running)
        scanning(mode)


def scanning(mode):
    global running
    global nominal

    global smr_real_posn_LTXYZ
    global benchinfo
    global label0, label1, label2

    if running == True:

        try:
            logger.info('Measuring reflector..')

            measured_LT = np.array(laser_tracker.measure())
            measurement = laser_tracker.measurement

            print('Measured SMR at LTXYZ: ', measured_LT)

            # Transform into the local optic coordinate system
            #
            measured_bench = benchinfo.transform_LT_to_bench(measured_LT)
            measured_zemax = benchinfo.transform_bench_to_Zemax(measured_bench)
            measured_optic = benchinfo.transform_Zemax_to_optic(measured_zemax)

            nominal_bench = benchinfo.transform_LT_to_bench(nominal)
            nominal_zemax = benchinfo.transform_bench_to_Zemax(nominal_bench)
            nominal_optic = benchinfo.transform_Zemax_to_optic(nominal_zemax)

            diffxyz = measured_optic - nominal_optic
            print("Measured position in bench coordsys: ", measured_bench)
            print("Measured position in zemax coordsys: ", measured_zemax)
            print("Measured position in optic coordsys: ", measured_optic)
            print("Position error in optic coordsys: ", diffxyz)

            ltlog(retro_id = -1, mode=mode, coordsys='LT',    xyz=measured_LT, err=[measurement.dStd1, measurement.dStd2, measurement.dStd3],redis_connection=myredis)
            ltlog(retro_id = -1, mode=mode, coordsys='bench', xyz=measured_bench,redis_connection=myredis)
            ltlog(retro_id = -1, mode=mode, coordsys='zemax', xyz=measured_zemax, zemax_file=zemax_prescrip_fname,zemax_surf=int(surface_number.get()),redis_connection=myredis)
            ltlog(retro_id = -1, mode=mode, coordsys='optic', xyz=measured_optic, zemax_file=zemax_prescrip_fname,zemax_surf=int(surface_number.get()),redis_connection=myredis)


            if mode == 'smr':
                label0.set( 'X: {:8.3f} mm | '.format(diffxyz[1]))
                label1.set( 'Y: {:8.3f} mm | '.format(diffxyz[0]))
                label2.set( 'Z: {:8.3f} mm'.format(diffxyz[2]))
                ltlog(retro_id = -1, mode=mode, coordsys='diff' , xyz=[diffxyz[1],diffxyz[0],diffxyz[2]],redis_connection=myredis)
            else:
                # Angular error of mirror
                hangle = np.degrees(diffxyz[1] / measured_optic[2] / 2) 
                vangle = np.degrees(diffxyz[0] / measured_optic[2] / 2)

                # Distance error at optic vertex
                dist   = (diffxyz[2]  + (measured_optic[0] - nominal_optic[0]) / nominal_optic[2] * nominal_optic[0] + (measured_optic[1] - nominal_optic[1]) / nominal_optic[2] * nominal_optic[1]) / 2
                label0.set( 'H: {:.5f} deg | '.format(hangle))
                label1.set( 'V: {:.5f} deg | '.format(vangle))
                label2.set( 'Z: {:.3f} mm'.format(dist))
                ltlog(retro_id = -1, mode=mode, coordsys='diff' , xyz=[hangle,vangle,dist],redis_connection=myredis)
        except:
            print('Error making measurement')

        tk.after(500, lambda: scanning(mode))


def stopMirrorAdjust():
    global running
    running = False

def initialize(laser_tracker):
    laser_tracker.initialize()
    reset()

running = False  # Global Flag

# permanent SMR measurements for SVD Rotation
#
permanent_smrs = np.loadtxt("permbench.txt")
if len(permanent_smrs) != len(config['smrs']):
    raise Exception("Configuration mismatch: %s contains %d SMRs but permbench.txt has %d lines" % ( config_file, len(config['smrs']), len(permanent_smrs)))

#
myredis = None
if config['redis']['use']:
    myredis = redis.Redis(host=config['redis']['host'])

# gui
tk = Tk()
tk.title("Mirror Alignment")

tk.columnconfigure(0, weight=1)
tk.rowconfigure(0, weight=1)

tk.option_add("*font", "lucida 12")

zemax_prescrip_fname = 'No Zemax File Selected!'
surface_number = StringVar()

# arrays for retro coordinates
cesapi_positions = np.loadtxt("tempbenchcoords.txt")
if len(cesapi_positions) != len(config['smrs']):
    raise Exception("Configuration mismatch: %s contains %d SMRs but tempbenchcoords.txt has %d lines" % ( config_file, len(config['smrs']), len(cesapi_positions)))

trackerpilot_positions = np.loadtxt("tempbenchcoords.txt")

nominal = 1

smr_real_posn_LTXYZ = 1
benchinfo = 1

laser_tracker = lt.LaserTracker(lt_host, logger)

# initialize button
Button(tk, text="Initialize", command=lambda: initialize(laser_tracker), bg="cyan").grid(
    column=0, row=0, sticky=(W, E))

# reconnect button
Button(tk, text="Reconnect", command=laser_tracker.reconnect, bg="cyan").grid(
    column=1, row=0, sticky=(W, E))

# smr widgets
smr_frame = LabelFrame(tk, text="MeasureSMRs")
smr_frame.grid(column=0, row=1, sticky=(N, W, E, S))
smr_frame.columnconfigure(0, weight=1)
smr_frame.rowconfigure(0, weight=1)

# Build a button for each of the SMRs, which are defined in the yaml config file.
smrs = config['smrs']
for smr_id in smrs.keys():
    smr = smrs[smr_id]
    measurecommand = partial(measureRetro, laser_tracker, index=smr_id)
    Button(smr_frame, text=smr['name'], command=measurecommand).grid(column=smr['layout_column'], row=smr['layout_row'], sticky=(W, E))


# coords widgets
coord_frame = LabelFrame(tk, text="Coordinate Options")
coord_frame.grid(column=0, row=2, sticky=(N, W, E, S))
coord_frame.columnconfigure(0, weight=1)
coord_frame.rowconfigure(0, weight=1)

Button(coord_frame, text="Load Coordinates", command=load).grid(column=0, row=0, sticky=(W, E))
Button(coord_frame, text="Save Coordinates", command=save).grid(column=1, row=0, sticky=(W, E))
Button(coord_frame, text="Reset Coordinates", command=reset).grid(column=2, row=0, sticky=(W, E))

# zemax widgets
zemax_frame = LabelFrame(tk, text="Zemax Options")
zemax_frame.grid(column=1, row=1, sticky=(N, W, E, S))
zemax_frame.columnconfigure(0, weight=1)
zemax_frame.rowconfigure(0, weight=1)

Button(zemax_frame, text="Select Zemax File", command=fileSelect).grid(column=0, row=0, sticky=(W, E))
Label(zemax_frame, text="Surface Number =").grid(column=1, row=0, sticky=(W, E))
Entry(zemax_frame, textvariable=surface_number).grid(column=2, row=0, sticky=(W, E))
Label(zemax_frame, text=(zemax_prescrip_fname), fg="Red").grid(columnspan=3, row=1, sticky=(W, E))

# adjust widgets
adjust_frame = LabelFrame(tk, text="Alignment and Positioning")
adjust_frame.grid(column=1, row=2, sticky=(N, W, E, S))
adjust_frame.columnconfigure(0, weight=1)
adjust_frame.rowconfigure(0, weight=1)

Button(adjust_frame, text="Run Mirror Adjust", command=lambda: Adjust("mirror"), bg="green").grid(column=0, row=0,
                                                                                                  sticky=(W, E))
Button(adjust_frame, text="Run SMR Adjust", command=lambda: Adjust("smr"), bg="green").grid(column=1, row=0,
                                                                                            sticky=(W, E))
Button(adjust_frame, text="STOP", command=stopMirrorAdjust, bg="red").grid(column=2, row=0, sticky=(W, E))

text_frame = Frame(tk)
text_frame.grid(columnspan=2, row=3)
text_frame.columnconfigure(0, weight=1)
text_frame.rowconfigure(0, weight=1)

# These labels will hold the position errors
#
label0 = StringVar(tk, value = '---')
label1 = StringVar(tk, value = '---')
label2 = StringVar(tk, value = '---')
Label(text_frame, textvariable = label0, font=("lucida", 48)).grid(column=0, row=0)
Label(text_frame, textvariable = label1, font=("lucida", 48)).grid(column=1, row=0)
Label(text_frame, textvariable = label2, font=("lucida", 48)).grid(column=2, row=0)

for child in smr_frame.winfo_children():
    child.grid_configure(padx=5, pady=5)
for child in coord_frame.winfo_children():
    child.grid_configure(padx=5, pady=5)
for child in zemax_frame.winfo_children():
    child.grid_configure(padx=5, pady=5)
for child in adjust_frame.winfo_children():
    child.grid_configure(padx=5, pady=5)

tk.mainloop()
