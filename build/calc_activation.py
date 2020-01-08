#!/usr/bin/env python
#fit activation energy and esitmate error from multiple calculations
#This fitting process works better for lowT; not so good on high T, as Ea is much more sensitive with same error of speed at high T

#Output
#$DIRECTION : x,y,z,r (r^2=x^2+y^2+z^2)
#Note diffusion coefficient is /6 for 3D case (r) or /2 for 1D case (x,y,z)

#$(Tempearture)-combine-$(DIRECITON).dat : time-distance at given temperature and direction, in time-distance(avg)-distance(min)-distance(max)-distance(error) format
#all-combine-$(DIRECTION).dat : time-distance(avg) for all temperatures, separated by empty lines
#speed-$(DIRECTION).dat : temperature-speed for all iterations
#speedloginv-$(DIRECTION).dat : temperature-speed for all iterations, in 1/temperature and ln(speed)
#ea-$(DIRECTION).dat : Fitted Ea using all temperatures by logK- 1/T, avg / min / max /error for all iterations. Note this value is not useful if temperature range is too large
#ea_t_$(DIRECTION).dat : Fitted Ea, using spline for logK - 1/T for all temperatures
#all-speed-avg.dat : temperature-speed , averaged for all iterations, for x/y/z/r




import os
import re
import sys
import platform
from math import *
import numpy as np
import shutil
from scipy.optimize import curve_fit
from scipy.interpolate import InterpolatedUnivariateSpline
from subprocess import check_call, Popen, PIPE
from argparse import ArgumentParser

parser = ArgumentParser(description="Calculate diffusion activation energy for given lattice")
parser.add_argument("case",nargs=1, help="Case name, a $CASE.in is required as input file")
parser.add_argument("-o", dest="dir_output", help="Output directory; by default the same as case name")
parser.add_argument("--exe",dest="st_exe", default=None, help="Executable path (abosulute path preferred), by default the 'randomwalk' in the same directoy")
parser.add_argument("--debug", dest="st_debug", default=None, required=False, help="The debugger to attach to each process")
parser.add_argument('-n',dest="n_repeat",type=int,default=10,required=False,help="Number of times runnning the whole simulation")
parser.add_argument('--np',dest="n_process",type=int,default=1,required=False,help="Number of parallel tasks")
args = parser.parse_args()

c0 = 1 #Prefactor
eV = 1.60217662e-19
kb = 1.38064852e-23 



st_exe = args.st_exe
if (st_exe is None):
    st_exe = os.path.join(os.path.dirname(os.path.realpath(__file__)), "randomwalk")

if (not os.path.exists(st_exe)):
    raise ValueError("Cannot find executable of randomwalk")

prefix = args.case[0]
dir_output = args.dir_output
if (dir_output is None):
    dir_output = prefix

print("Prefix : %s" % prefix)
filename_in = "%s.in" % prefix
filename_latt = "%s.latt" % prefix
re0 = re.compile("%s-t-(.+).dat" % prefix)

n_process = args.n_process #Number of processors
n_repeat = args.n_repeat #Repeat calculations
n_sample = 100
n_trajectory = 10000

if (os.path.exists(dir_output)):
    print("Folder %s exists, directly read" % dir_output)
else:
    os.mkdir(dir_output)
#Backup input files
    dir0 = os.getcwd()
    for filename in (filename_in, filename_latt):
        shutil.copy(os.path.join(dir0, filename), os.path.join(dir_output, filename))

    n_batch =  n_repeat / n_process
    for ix_iter in range(n_batch):
        list_proc = []
        list_ix_repeat = [x+1 for x in range(ix_iter * n_process, min(n_repeat, (ix_iter+1) * n_process))]
        for ix_repeat in list_ix_repeat:
            dir1 = os.path.join(dir_output, "%i" % ix_repeat)
            print("Iteration %i" % ix_repeat)
            os.mkdir(dir1)
            for filename in (filename_in, filename_latt):
                shutil.copy(os.path.join(dir0, filename), os.path.join(dir1, filename))
            os.chdir(dir1)
            if (platform.system() == "Windows"):
                proc = Popen(["../../../x64/Release/RandomWalk.exe",prefix], stderr=PIPE)
            else:
                if (args.st_debug is not None):
                    proc = Popen([args.st_debug, "--args", st_exe, prefix], stdin=PIPE, stdout=PIPE, stderr=PIPE,
                            close_fds=True)
                    proc.stdin.write("run\nbt\nquit\ny\n")
                    proc.stdin.close()
                else:
                    proc = Popen([st_exe, prefix], stdout=PIPE, stderr=PIPE)
            os.chdir(dir0)
            list_proc.append(proc)
        for proc, ix_repeat in zip(list_proc, list_ix_repeat):
            exitcode = proc.wait()
            if (exitcode != 0):
                print("Error!")
                print(proc.stdout.read())
                print(proc.stderr.read())
                raise ValueError("Non-zero exitcode for %i" % ix_repeat)
            elif (args.st_debug is not None):
                print("Process %i" % ix_repeat)
                print(proc.stdout.read())
                print(proc.stderr.read())


list_data_repeat = []
for ix_repeat in range(n_repeat):
    dir1 = os.path.join(dir_output, "%i" % (ix_repeat+1))
    list_data = []
    for dirroot, dirnames, filenames in os.walk(dir1):
        for filename in filenames:
            m = re0.match(filename)
            if (m is not None):
                list_data.append([float(m.group(1)), np.loadtxt(os.path.join(dirroot, filename))])
        break
    list_data.sort(key=lambda x:x[0])
    list_temperature = np.asarray([x[0] for x in list_data])
    list_data_repeat.append([x[1] for x in list_data])

#Check consistency
ar_len = [len(x) for x in list_data_repeat]
if (len(set(ar_len)) != 1):
    print(ar_len)
    raise ValueError("Temperature not consistent between repeat runs")
ar_data_xyz = np.asarray(list_data_repeat)

#Subscript: repeat, temperature, timestep, time/x**2/y**2/z**2 
#Reform to temperature, timestep, xyzr, time/avg/min/max (min max are min/max in all repeats)
ar_data_r = ar_data_xyz[:,:,:,1] + ar_data_xyz[:,:,:,2] + ar_data_xyz[:,:,:,3]
ar_data_xyzr = np.append(ar_data_xyz, ar_data_r[:,:,:,np.newaxis], axis=3)

list_speed_alldirection = [list_temperature]

for ixyzr, xyzr in [
        (1, "x"),
        (2, "y"),
        (3, "z"),
        (4, "r"),
        ]:
    print("Direction : %s" % xyzr)

    ar_data = ar_data_xyzr[:,:,:,[0,ixyzr]]

    ar_avg = ar_data[:,:,:,1].mean(0)
    ar_min = ar_data[:,:,:,1].min(0)
    ar_max = ar_data[:,:,:,1].max(0)
#First step is always zero ; manually 
    ar_err = np.maximum(ar_max - ar_avg, ar_avg - ar_min)[:,1:] / ar_avg[:,1:]
    ar_err = np.concatenate([ar_err, np.zeros((ar_err.shape[0],1))], axis=1)

    ar_data_repeat = ar_data
    ar_data = np.stack((ar_data[0,:,:,0], ar_avg, ar_min, ar_max, ar_err), axis=-1)
#print(ar_data.shape)

#Write all combined data
    dic_temperature_error = {}
    for i1, temperature in enumerate(list_temperature):
        np.savetxt(os.path.join(dir_output, "%s-combine-%s.dat" % (temperature,xyzr)), ar_data[i1,:,:])
        #Estimate error; only use first 75%
        nstep = ar_data.shape[1]
        sample = ar_data[i1,nstep*1/4:nstep*3/4,1:]
        error = np.mean(np.maximum(sample[:,0]-sample[:,1], sample[:,2]-sample[:,0])/sample[:,0])
#   print("Error for temperature %.3f : %.2f%%" % (temperature, 
#       100*error))
        dic_temperature_error[temperature] = error

#t-R^2 averaged, for all T
    with open(os.path.join(dir_output, "all-combine-%s.dat" % xyzr), 'w') as f:
        for i1, temperature in enumerate(list_temperature):
            np.savetxt(f, ar_data[i1, :, 0:2])
            f.write("\n")

    error_max = max(dic_temperature_error.values())
#print("Maximum error : %.2f%%" % (100 * error_max))

#Convert trajectories to activation energy
#Get speed constant : <L^2> / t

    def f_t_L2(t, c):
        return c*t

    def f_T_c(T, Ea, A):
        return A*np.exp(-Ea * (eV / kb / T))

    def f_T_logc(T, Ea, A):
        return A + -Ea * (eV / kb / T)


#Delete a percentage as warm up (not included in the fitting
    ratio_warmup = 0.25
    n_step = int(ar_data_repeat.shape[2] * ratio_warmup)

    list_speed_iter = []

    for i in range(ar_data_repeat.shape[0]):
        list_speed = []
        for i1, temperature in enumerate(list_temperature):
            data = ar_data_repeat[i, i1, :, :]
            list_speed.append([temperature, curve_fit(f_t_L2, data[n_step:,0], data[n_step:,1], p0=[1.0])[0][0]])

        list_speed.sort(key=lambda x:x[0])
        ar_speed = np.asarray(list_speed)
        #for x1, x2 in ar_speed:
        #    print("%f %e" % (x1, x2))
        list_speed_iter.append(ar_speed)

#Estimate speed error
#Subscript : iter, temperature, 2 (first temperature value and second speed)
    ar_speed_iter = np.asarray(list_speed_iter)[:,:,1]
    ar_speed_mean = ar_speed_iter.mean(axis=0)
    ar_speed_error1 = (ar_speed_iter.max(axis=0) - ar_speed_mean)/ar_speed_mean
    ar_speed_error2 = -(ar_speed_iter.min(axis=0) - ar_speed_mean)/ar_speed_mean
    ar_speed_error = np.stack([ar_speed_error1, ar_speed_error2]).max(axis=0)
    print("Speed estimated error:")
    for i, temperature in enumerate(list_temperature):
        print("Temperature %7.2f K :  Speed Error %.2f%%" % (temperature, ar_speed_error[i]*100))
    list_speed_alldirection.append(ar_speed_mean)

#Join speed two columns: (temperature , speed)
    ar_speed_all = np.stack([list_speed_iter[0][:,0]] + [x[:,1] for x in list_speed_iter], axis=1)
    np.savetxt(os.path.join(dir_output, "speed-%s.dat" % xyzr), ar_speed_all)

    ar_speed_all[:,0] = 1.0/ ar_speed_all[:,0]
    ar_speed_all[:,1] = np.log(ar_speed_all[:,1])

    np.savetxt(os.path.join(dir_output, "speedloginv-%s.dat" % xyzr), ar_speed_all)

#Fit Ea
    list_ea = []
    for ar_speed in list_speed_iter:
        popt, pcov =  curve_fit(f_T_logc, ar_speed[:,0], np.log(ar_speed[:,1]), p0=[0.2, 1])
        Ea = popt[0]
        list_ea.append(Ea)

        
    ar_ea = np.asarray(list_ea)
    ea_stat = np.asarray([ar_ea.mean(), ar_ea.min(), ar_ea.max(),
            -(ar_ea.mean() - ar_ea.min())/ar_ea.mean()/10,
            -(ar_ea.mean() - ar_ea.max())/ar_ea.mean()/10, 
            ])
    print("Activation energy : %.2f meV (Range %.2f~%.2f %.2f%%~%.2f%%)" % tuple(ea_stat*1000))
    np.savetxt(os.path.join(dir_output, "ea-%s.dat" % xyzr), ea_stat)

#Fit ea at different temperature range
    fit_order = 3 #Need high order to converge
    list_ea = []
    for ar_speed in list_speed_iter:
        invT = 1.0/ar_speed[:,0]
        logC = np.log(ar_speed[:,1])
#Note : Keep increasing x values in fitting!
#Polynomial fitting (not good)
#   p = np.polyder(np.polyfit(invT[::-1], logC[::-1], fit_order))
#   ar_ea = np.polyval(p, invT)
#Spline fitting
        p = InterpolatedUnivariateSpline(invT[::-1], logC[::-1], k=fit_order)
        p = p.derivative(n=1)
        ar_ea = [p(x) for x in invT]

        list_ea.append(ar_ea)
    ar_ea = np.asarray(list_ea)

    def get_stat(ar):
        amean = ar.mean(axis=0)
        amax = ar.max(axis=0)
        amin = ar.min(axis=0)
        return np.asarray([amean, np.maximum(np.abs(amax-amean), np.abs(amin-amean))/amean]).T

    ea_stat = get_stat(ar_ea)

    print("Temperature dependent Ea by %i-degree fit" % fit_order)
    for temperature, row in zip(list_temperature, ea_stat):
        mean = row[0]
        err = row[1]
        print("Temperature %7.2f K : %.2f meV (Error %.2f%%)" % (temperature, mean * kb / eV * 1000,   err*100))
    np.savetxt(os.path.join(dir_output, "ea_t_%s.dat" % xyzr), np.concatenate((list_temperature[:,np.newaxis], ea_stat), axis=1))

np.savetxt(os.path.join(dir_output, "all-speed-avg.dat"), np.asarray(list_speed_alldirection).T)
