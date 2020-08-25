//gcc -std=c99 -o go main2.c PhaseCalc.c -lm

#ifndef __hirax_phase_calc
#define __hirax_phase_calc

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#ifndef C
#define C 300000000.0
#endif

/* Convert from decimal degrees to radians */
#define degrad(x)    ((x)*M_PI/180.)

/* Convert from radians to decimal degrees */
#define raddeg(x)    ((x)*180./M_PI)

//Single precision
//****************************************************************************************************//
typedef struct {
    
    float real; //real part of the compex number
    float imag; //imaginary part of the compex number
}complex;


//****************************************************************************************************//
/**
 * @struct ZenAzi - structure containing angles zenith and azimuth angles
 *za -zenith-angle in radians
 *az - azimuth angle in radians, measured from East(x) to North(y).
 **/
typedef struct {
    float alt, az;
}AltAz_cordinates;


//****************************************************************************************************//
/**
 * @struct  - structure containing angles hour angle an declination
 *ha - hour angle in radians
 *dec - declination in radians
 **/
typedef struct {
    float ha, dec;
}HADEC_cordinates;



//****************************************************************************************************//
/**
 *@struct XYZ - structure containing direction cosines coordinates (x, y, z)
 * x points (East) direction, in metres
 * y points (North) direction, in metres
 * z points (Up), in metres
 **/

typedef struct {
    float x, y, z;
}XYZ_cordinates;

//****************************************************************************************************//
/**
 * @brief - that converts spherical coordinates to enu direction cosines
 *
 ***In parameters***
 *-------------------
 *  theta - polar (zenith) angle in radians
 *  phi  - co-azimuth angle in radians
 *   (phi = 0 ==> x, y = 1,0)
 *   (phi = 90 ==>x,y = 0,1
 
 ***Out Parameters***
 *----------------
 * structure containing direction cosines coordinates (x, y, z) relative to zenith.
 **/

XYZ_cordinates *convert_azel_to_dircos(AltAz_cordinates thtph);



//****************************************************************************************************//
/**
 * @func that converts hour angle, declination to enu direction cosines
 *
 *@In parameters
 *----------------
 *  ha - hour angle in radians.
 *  dec - declination angle  in radians.
 *  latitude - observers geodetic latitude in radians
 
 *@Out Parameters
 *----------------
 * structure containing direction cosines coordinates in horizontal system (Xh, Yh, Zh).
 xh - the east-west direction horizontal cordinates.
 yh - the north-south direction horizontal cordinates.
 zh - the up direction horizontal cordinates.
 **/
XYZ_cordinates *convert_hadec_to_enu(HADEC_cordinates hadec, float lat);


//****************************************************************************************************//
/**
 *@brief A function that convert Horizontal coordinates to spherical cordinates (az, alt).
 *@param[In]: struct of horizontal coordinates (x,y,z)
 
 *@param[out]: struct spherical coordinates (az, alt).
 az - azimuth angels in radians.
 alt - altitude angels in radians.
 **/
AltAz_cordinates *convert_enu_to_az_alt(XYZ_cordinates *enu_cordinates);

//****************************************************************************************************//



void compute_antenna_phases(TABs_inform *beams_param, float *ant_gains, complex **OutPhase_east, complex **OutPhase_north , float lst, float freq);

#endif


