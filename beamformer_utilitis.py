import numpy as np
from unitsconvert import *
from astropy import constants as cons
c = cons.c.value
from Coordinates import *

def get_XY_steering_vectors(xPos, yPos, freq, azimuth, altitude):
    
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
    [X,Y,Z] = altaz2dircos(theta_grid, phi_grid)

    psi_x = X[np.newaxis,:,:] * xPos[:, np.newaxis, np.newaxis]
    psi_y = Y[np.newaxis,:,:] * yPos[:, np.newaxis, np.newaxis]

    xSteering_Vectors = np.exp(-1j * phaseConstant * psi_x)
    ySteering_Vectors = np.exp(-1j * phaseConstant * psi_y)
    
    return [xSteering_Vectors, ySteering_Vectors]



def TABs_steering_vectors(xPos, yPos, freq, beam_az, beam_alt):
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
        """
    lamb = freq/c
    phaseConstant = (2 * math.pi)/lamb
    [l_source,m_source,n_source] = altaz2dircos(beam_az, beam_alt)
    tmp_x = l_source[:,np.newaxis] * xPos
    tmp_y = m_source[:,np.newaxis] * yPos
    xbeam_Steering_Vectors = np.exp(1j * phaseConstant * tmp_x)
    ybeam_Steering_Vectors = np.exp(1j * phaseConstant * tmp_y)
    return [xSource_Steering_Vectors, ySource_Steering_Vectors]
