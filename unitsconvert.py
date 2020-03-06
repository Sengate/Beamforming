import numpy as np
import math
import matplotlib.pyplot as plt

#================================================================================================
#                               Units conversions
#=----------------------------------------------------------------------------------------------------

PI = math.pi

def arcsec2deg(angle):
    """convert arc-seconds to degrees"""
    return angle/3600.0

def arcmin2deg(angle):
    """convert arcminutes to degree"""
    return angle/60

def DegMinSec2rad(angle):
    "Convert angle given in deg + arcmin + arcsec to degree fraction"
    deg, arcmin, arcsec = angle
    degree = deg +  arcmin2deg(arcmin) + arcsec2deg(arcsec)
    return deg2rad(degree)

def DegMinSecDeg(angle):
    "Convert angle given in deg + arcmin + arcsec to degree fraction"
    deg, arcmin, arcsec = angle
    degree = deg +  arcmin2deg(arcmin) + arcsec2deg(arcsec)
    return degree

def deg2rad(deg):
    """convert degrees to radians"""
    return deg * (PI/180.0)

def rad2deg(rad):
    """convert radians to degrees"""
    return rad * (180.0/PI)

def time2dec(Time):
    """convert time in hours, minutes, and seconds to decimal format"""
    hours, minutes , seconds = Time
    return hours + (minutes/60.0)+ (seconds/(3600.0))

def hr2rad(Time_dec):
    """Convert time hours decimal format to radians
        
        :1rad = 15/(180/pi)
        """
    return (Time_dec * 15.) * PI/180.0


def rad2hr(rad):
    """convert radians to hours"""
    return rad * 12.0/PI


def hr2deg(hr):
    """Convert time in hours decimal format to degress
        
        :1hr = 15 degrees """
    return hr * 15.0

def deg2hr(deg):
    return deg/15.0

