#include <stdio.h>
#include <stdlib.h>

#include "PhaseCalc.h"

XYZ_cordinates *get_XYZ(float x,
             float y,
             float z
             )
{
    XYZ_cordinates *pr = (XYZ_cordinates*)malloc(sizeof(XYZ_cordinates));
    pr->x = x;
    pr->y = y;
    pr->z = z;
    
    return pr;
    
    free(pr);
}

/*--------------------------------------------------------------------------*/
//Convert the theta, phi coordinate to direction cosines
XYZ_cordinates *convert_azel_to_dircos(AltAz_cordinates thtph){
    
    float theta = thtph.alt;
    float phi = thtph.az;
    float x = sinf(theta) * cosf(phi);
    float y = sinf(theta) * sinf(phi);
    float z = cosf(theta);
    
    return get_XYZ(x, y, z);
}


/*--------------------------------------------------------------------------*/
//Convert hadec coordinates to ENU given the telescope latitude and beams (ha, dec)
XYZ_cordinates *convert_hadec_to_enu(HADEC_cordinates hadec, float obs_latitude){
    float Xh,Yh, Zh;
    float ha = hadec.ha;
    float dec = hadec.dec;
    Xh = -sinf(dec) * cosf(obs_latitude) + sinf(obs_latitude) * cosf(dec) * cosf(ha);
    Yh = sinf(ha) * cosf(dec);
   // Yh = sinf(dec) * cosf(obs_latitude) - sinf(obs_latitude) * cosf(dec) * cosf(ha);
    Zh = cosf(obs_latitude) * cosf(dec) * cosf(ha) + sinf(obs_latitude) * sinf(dec);
    
    return get_XYZ(Xh, Yh, Zh);
}


/*---------------------------------------------------------------------------------*/
//convert the ENU coordinates to altaz
AltAz_cordinates *convert_enu_to_az_alt(XYZ_cordinates *enu_cordinates){
    
    float az = atan2f(enu_cordinates->y, enu_cordinates->x); //azimuth angle
    float xx = sqrtf(enu_cordinates->x * enu_cordinates->x + enu_cordinates->y*enu_cordinates->y);
    float el = atan2f(enu_cordinates->z, xx); //elevation
   //float el = asinf(enu_cordinates->z);
    
    AltAz_cordinates *altaz = (AltAz_cordinates*)malloc(sizeof(AltAz_cordinates));
    altaz->az = az;
    altaz->alt = el;
                     
    return altaz;
    
    free(altaz);
}


/*****************************************BEams Phases Calculations*********/

// beams information
typedef struct{
    float *RAs;     //beams right acessions
    float *DECs;     // beams declinations
}TABs_inform;


void compute_antenna_phases(TABs_inform *beams_param, float *ant_gains, complex **OutPhase_east, complex **OutPhase_north , float lst, float freq){
    //lst in radians
    //beams positions in radians
    
   
    float *ha = (float*)malloc(sizeof(float) * nbeams);
    float *beams_azimuths = (float*)malloc(sizeof(float) * nbeams);
    float *beams_altitudes = (float*)malloc(sizeof(float) * nbeams);
    
    float lam = C/freq;
    
    for (int ibeam=0; ibeam < nbeams; ibeam++){ //loop over beams
        
        //get hour angles from RAs.
        ha[ibeam] = lst - beams_param->RAs[ibeam];
        HADEC_cordinates EQU_cord = {ha[ibeam], beams_param->DECs[ibeam]};
        
        //convert hadec to enu coordinates (x, y, z).
        XYZ_cordinates *ENU_cord = convert_hadec_to_enu(EQU_cord, obs_latitude);
        
        //get beam azimuths and altitude angles
        AltAz_cordinates *azalt = convert_enu_to_az_alt(ENU_cord);
        beams_azimuths[ibeam] = azalt ->az;
        beams_altitudes[ibeam] = azalt ->alt;
        AltAz_cordinates AzAlt = {beams_azimuths[ibeam], beams_altitudes[ibeam]};
        
        //convert altaz coordinates to direction cosines
        XYZ_cordinates *DirecCos = convert_azel_to_dircos(AzAlt);
        
        for (unsigned int iant = 0; iant < nants; iant++){ //loop over antennas

            //East_west pojection phase vector
            float angle_EW = 2 * M_PI/lam * (iant * east_antSpacing) * DirecCos->x;
            OutPhase_east[ibeam][iant].real =cos(angle_EW) * ant_gains[iant];
            OutPhase_east[ibeam][iant].imag =-sin(angle_EW) * ant_gains[iant];
            
            //North_south projection phase vector
            float angle_NS = 2 * M_PI/lam * (iant * north_antSpacing) * DirecCos->y;
            OutPhase_north[ibeam][iant].real =cos(angle_NS) * ant_gains[iant];
            OutPhase_north[ibeam][iant].imag =-sin(angle_NS) * sant_gains[iant];
            
            
        }
    }
}
    
    
    
    

