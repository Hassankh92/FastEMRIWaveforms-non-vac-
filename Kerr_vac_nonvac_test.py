# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 00:12:59 2022

@author: hassa
"""
import numpy as np
import scipy as sp
import scipy.special
import numpy as np
import math
import matplotlib.pyplot as plt


import warnings

from matplotlib import pyplot as plt
from few.trajectory.inspiral import EMRIInspiral
from few.utils.baseclasses import SchwarzschildEccentric, Pn5AAK, ParallelModuleBase, KerrCircular
from few.waveform import FastSchwarzschildEccentricFlux, AAKWaveformBase, Pn5AAKWaveform, SlowSchwarzschildEccentricFlux, KerrCircularFlux, MigrationTorqueKerrCircularFlux, ScalarCloudKerrCircularFlux

from few.summation.aakwave import AAKSummation
from few.utils.utility import get_mismatch
from few.waveform import   GenerateEMRIWaveform
from few.summation.interpolatedmodesum import CubicSplineInterpolant
from few.utils.constants import MTSUN_SI, YRSID_SI, Pi
from few.waveform import SchwarzschildEccentricWaveformBase, RelativisticKerrCircularWaveformBase
from few.utils.baseclasses import SchwarzschildEccentric
from few.trajectory.inspiral import EMRIInspiral
from few.amplitude.interp2dcubicspline import Interp2DAmplitude, Interp2DAmplitudeKerrCircular
from few.summation.interpolatedmodesum import InterpolatedModeSum, CubicSplineInterpolant, InterpolatedModeSumKerrCircular
from few.summation.directmodesum import DirectModeSum
from few.utils.ylm import GetYlms
from few.amplitude.romannet import RomanAmplitude

try:
    import cupy as xp

    gpu_available = True

except (ModuleNotFoundError, ImportError) as e:
    import numpy as xp

    warnings.warn(
        "CuPy is not installed or a gpu is not available. If trying to run on a gpu, please install CuPy."
    )
    gpu_available = False


traj_few = EMRIInspiral(func="SchwarzEccFlux")
print("old few traj ran",traj_few,"\n")

traj_Kerr = EMRIInspiral(func="KerrCircFlux")
print("Kerr traj ran",traj_Kerr,"\n")


traj_migration = EMRIInspiral(func="MigTorqKerrCircFlux")
traj_cloud = EMRIInspiral(func="CloudKerrCircFlux")


T = 4.0 # years
dt = 10.0
# set initial parameters
M = 1e6
mu = 5e1
p0 = 12.0
e0 = 0.0
a0 = 0.5
x0 = 1.0
Y0 = 1.0


# sample values for migration params
A = 1e-5
nr = 8

# a sample value for the cloud
alpha =  0.16
Mb_M = 0.05



Kerr_result = traj_Kerr(M,mu,a0,p0,e0,Y0,T=T, dt=dt)
migration_result = traj_migration(M,mu,a0,p0,e0,Y0,A, nr,T=T, dt=dt)
cloud_result = traj_cloud(M,mu,a0,p0,e0,Y0,alpha,Mb_M,T=T, dt=dt)

t_Kerr  = Kerr_result[0]
t_migration = migration_result[0]
t_cloud = cloud_result[0]
print("Kerr final time: ", t_Kerr[-1], "Migration final time: ", t_migration[-1], "Cloud final time: ", t_cloud[-1])


p_Kerr = Kerr_result[1]
p_migration = migration_result[1]
p_cloud = cloud_result[1]
print("Kerr final p: ", p_Kerr[-1], "Migration final p: ", p_migration[-1], "Cloud final p: ", p_cloud[-1])
