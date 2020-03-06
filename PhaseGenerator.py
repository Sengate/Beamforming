from astropy import units
from astropy import constants as cons
import numpy as np
import matplotlib.pyplot as plt
import math
import astropy.units as u
from unitsconvert import *

c = cons.c.value

def AltAz_to_Dircos(azimuth, altitude, from_zenith=False):
    
    """
        Compute direction cosines (X,Y, Z) from azimuth and altitude
        
        Parameters
        ----------
        
        azimuth: the azimuth angle measured clockwise for North-East wards [rad]
        Can be a tuple, list or np.array
        
        altitude: the altitude
        Can be a tuple, list or np.array
        
        return
        -------
        
        X - pointed towards north-south point
        Y - chosen such that the system is left-handed (points East-West)
        Y - points towards the zenith
        
        """

    if from_zenith:
        altitude = math.pi/2 - altitude
    #else:
    #   altitude  =  altitude
    X = np.cos(azimuth) * np.cos(altitude)
    Y = np.sin(azimuth) * np.cos(altitude)
    Z = np.sin(altitude)
                
    return [X,Y,Z]


def hadec_to_HorizDircos(ha, dec, lat):
    
    """
        Compute the Direction cosine of the object in from the hour angle and declination
        
        This function gives the relationship between the  equatorial coordinates (Xe,Ye,Ze) and
        altazimuth coordinates (Xh,Yh,Zh)
        
        prameters
        ----------
        ha: Hour Angle [rad]
        Can be tuple, list or np.array
        
        dec: Declination [rad]
        Can be tuple, list or np.array
        
        lat: Obsevers latitude [rad]
        float
        
        returns:
        --------------
        List of the Horizontal coordinates (Xh, Yh, Zh)
        
        """
    
    Xh = np.sin(dec) * np.cos(lat) - np.sin(lat) * np.cos(dec) * np.cos(ha)
    Yh = -np.sin(ha) * np.cos(dec)
    Zh = np.cos(lat) * np.cos(dec) * np.cos(ha) + np.sin(lat) * np.sin(dec)
    return [Xh, Yh, Zh]


def Compute_AltitudeAzimuth(dircos):
    
    """ Compute the Altitude and Azimuth of objects, given its Direction cosines
        
        Parameters:
        dircos: List, Array or tuple of (Xh, Yh, Zh) [rad]
        
        Returns
        ---------
        alt: Altitude [rad]
        az: Azimuth [rad]
        
        """
    if not any([isinstance(dircos,list), isinstance(dircos, np.ndarray), isinstance(dircos, tuple)]):
        raise TypeError('thph should be .....')
    #if isinstance(dircos, tuple):
    
    dircos= np.array(dircos)
    
    X = dircos[0].ravel()
    Y = dircos[1].ravel()
    Z = dircos[2].ravel()
    
    alt = np.arcsin(Z)
    #az = np.arctan(Y/X)
    az = np.arctan2(Y, X)
    
    az[np.where(az<0)]+=2*np.pi
    
    alt = alt.reshape(dircos[0].shape)
    az = az.reshape(dircos[0].shape)
    
    return [alt, az]


def Get_HA(ra, lst):
    """compute the hour angle given the local sidereal time (LST) and Right Ascension(RA)
        
        Prarameters
        --------------------------------------
        _LST: local sidereal time [Hours]
        RA : right ascension [Hours]
        
        returns:
        (ha[hrs], ha[rad], ha[deg])
        """
    ra = rad2hr(ra)
    ha = lst - ra # HA in hours fraction
    ha_deg = np.mod(hr2deg(ha), 360.0) #take the modulus angle
    ha_RAD = deg2rad(ha_deg)
    #return HA in (hours, radians, degree)
    return ha_RAD


def Compute_XYsteering_Vectors(xPos, yPos, freq, azimuth, altitude):
    
    """
        Parameters:
        xPos        1xN of x-axis antenna positions:[m]
        yPos        1xN of y-axis antenna positions:[m]
        freq        Observing freq [Hz]
        azimuth     Azimuth angle measured in the XY-plane [rad]
        altitude    Altitude angles normal to the incident angle [rad]
        
        Output:
        xSteering_Vectors - A matrix of steering vectors in the x-direction with shape (az.shape, alt.shape, N)
        ySteering_Vectors - A matrix of steering vectors in the y-direction with shape (az.shape, alt.shape, N)
        """
    
    lamb = freq/c
    phaseConstant = (2 * math.pi)/lamb
    
    #get theta and phi grids
    theta_grid = altitude.reshape(-1,1)
    phi_grid = azimuth.reshape(1,-1)
    
    #get directional cosines
    [X,Y,Z] = AltAz_to_Dircos(theta_grid, phi_grid)
    psi_x = X[np.newaxis,:,:] * xPos[:,np.newaxis, np.newaxis]
    psi_y = Y[np.newaxis,:,:] * yPos[:,np.newaxis, np.newaxis]
    xSteering_Vectors = np.exp(-1j * phaseConstant * psi_x)
    ySteering_Vectors = np.exp(-1j * phaseConstant * psi_y)
    
    return [xSteering_Vectors, ySteering_Vectors]



def Compute_beams_SteeringVectors(xPos, yPos, freq, beam_az, beam_alt):
    """
        Parameters:
        xPos        1xN of x-axis antenna positions:[m]
        yPos        1xN of y-axis antenna positions:[m]
        freq        Observing freq [Hz]
        beam_az     TAB azimuth angle [rad]
        beam_alt    TAB altitude angle [rad]
        
        Output:
        xbeam_Steering_Vectors -
        ybeam_Steering_Vectors -
        
        shape(nbeams, nant)
        """
    lamb = freq/c
    phaseConstant = (2 * math.pi)/lamb
    [X_beam,Y_beam,Z_beam] = AltAz_to_Dircos(beam_az, beam_alt)
    tmp_x = X_beam[:,np.newaxis] * xPos
    tmp_y = Z_beam[:,np.newaxis] * yPos
    xbeam_Steering_Vectors = np.exp(1j * phaseConstant * tmp_x)
    ybeam_Steering_Vectors = np.exp(1j * phaseConstant * tmp_y)
    
    return [xbeam_Steering_Vectors, ybeam_Steering_Vectors]



def main(frequencies, xPos, yPos, beams_RAs, beams_Decs, LST_list=None, LSTmin=None,
         LSTmax=None,THETA=None,PHI=None,
         field_of_view=None, Longitude=None, Latitude=None, Telescope=None):
    
    
    """
        Parameters
        -----------
        frequencies - is observing frequencies in Hz.
        It can be list or a numpy array
        
        Ant_positions - a file containing (X,Y,Z) or (E,N,U) antenna positions
        
        beams_list - right ascessions and declinations of then beams (RAs, DECs) [radians]
        
        PHI - is the boresite beam azimuth angles in radians
        Should be a grid of azimuth angles
        
        THETA - is the boresite beam altitude angles in radians
        hould be a grid of altitude angles
        
        field_of_view - ((tht_lo, ph_lo),(tht_hi, phi_hi))
        
        LST_list - list of local sidereal times in fraction hours
        
        Longitude - array longitude in radians
        
        Latitude - array latitude in radians
        
        Telescope - telescope position (latitute, longitude, height)
        
        returns
        
        ZENITH_PHASEVECTORS - A dictionary containing phase vector for NS and EW direction, for each frequency
        and antenna
        
        BEAMS_PHASEVECTORS - A dictionary containing phase vectors for each beam, each frequecny, LST time step and antenna
        -------------
        
        """
    
    #get antenna parameters
    #xPos = Ant_positions[:,0]
    #yPos = Ant_positions[:,1]
    #zPos = Ant_positions[:,2]
    nant = len(xPos) #number of antennas
    
    #get beams parameters
    #beams_RAs = beams_list[:,0]
    #beams_Decs = beams_list[:,1]
    if len (beams_RAs) != len(beams_Decs):
        raise ValueError(' RAs and DECs must have equal length().')
    nbeams = len(beams_RAs)
    
    
    #get LST parameters
    if LST_list:
        lst_list = LST_list
    else:
        lst_list = []
        count = int(math.ceil((LSTmax - LSTmin) / delta))
        for i in range(0, count+1):
            lst_list.append(LSTmin + i*delta)
    ntime = len(lst_list)
    
    #get telescope parameters
    if Telescope:
        latitude, longitude, hieght = telescope
    else:
        latitude = Latitude
        longitude = Longitude
    
    
    
    #------------------------------------------------------------------------------------------#
    ##              Get Phases
    #------------------------------------------------------------------------------------------#
    ZenithSteeringVectors = []
    #has = []
    #ntime_BeamsAltitudes = []
    #ntime_BeamsAzimuths = []
    beams_steeringVectors = []


    #compute the zenith steering vectors for each frequency
    #get North-South & East-West Steering vectors matrix with shape = (nfreq,nants,az.shape, alt.shape)
    for ifreq, freq in enumerate(frequencies):
        ZenithSteeringVectors.append(Compute_XYsteering_Vectors(xPos, yPos, freq, PHI, THETA))
    ZenithSteeringVectors = np.asarray(ZenithSteeringVectors)
    
        #shape = (nfreq, nant, az.shape, alt.shape)
    NS_ZenithSteeringVectors = ZenithSteeringVectors[:,0,:,:,:]
    EW_ZenithSteeringVectors = ZenithSteeringVectors[:,1,:,:,:]
    
    
    #compute beams altitudes and azimuths for each time step and frequency
    #compute beams steering vectors for each time step and frequency
    #get North-South & East-West Steering vectors matrix with shape = (nfreq,ntimes,nbeams,nants)
    for ifreq, freq in enumerate(frequencies):
        
        beams_steeringVectors.append([])
        
        for itime, time in enumerate(lst_list):
            
            beams_RAs = np.asarray(beams_RAs)
            HAs = Get_HA(beams_RAs, time)
            #has.append(HAs)
            
            #get (Xh,Yh,Zh) from the hour angle and declination
            dircos = hadec_to_HorizDircos(HAs, beams_Decs, latitude)
            
            #get altaz from (Xh,Yh,Zh)
            beams_altaz = Compute_AltitudeAzimuth(dircos)
            
            beams_altitudes = beams_altaz[0]
            beams_azimuths = beams_altaz[1]
            
            #get altitudes and azimuths
            #ntime_BeamsAltitudes.append(beams_altitudes)
            #ntime_BeamsAzimuths.append(beams_azimuths)
            
            #get beam steering vectors
            beams_Stv=Compute_beams_SteeringVectors(xPos, yPos, freq, beams_azimuths, beams_altitudes)
            beams_steeringVectors[ifreq].append(beams_Stv)

    #ntime_BeamsAltitudes = np.asarray(ntime_BeamsAltitudes)
    #ntime_BeamsAzimuths = np.asarray(ntime_BeamsAzimuths)

    beams_steeringVectors = np.asarray(beams_steeringVectors)

    NS_beams_SteeringVectors = beams_steeringVectors[:,:,0,:,:]
    EW_beams_SteeringVectors = beams_steeringVectors[:,:,1,:,:]
    
    
    
    BEAMS_PHASEVECTORS = {}
    BEAMS_PHASEVECTORS['frequencies'] = frequencies
    BEAMS_PHASEVECTORS['times'] = lst_list
    BEAMS_PHASEVECTORS['NS_beams_SteeringVector'] = NS_beams_SteeringVectors
    BEAMS_PHASEVECTORS['EW_beams_SteeringVectors'] = EW_beams_SteeringVectors
    
    ZENITH_PHASEVECTORS = {}
    ZENITH_PHASEVECTORS['frequencies'] = frequencies
    ZENITH_PHASEVECTORS['NS_ZenithSteeringVectors'] = NS_ZenithSteeringVectors
    ZENITH_PHASEVECTORS['EW_ZenithSteeringVectors'] = EW_ZenithSteeringVectors
    
    
    return [ZENITH_PHASEVECTORS, BEAMS_PHASEVECTORS]


