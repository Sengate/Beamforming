
#ifndef hirax_beamform_hadec_to_enu
#define hirax_beamform_hadec_to_enu


#include <stdlib.h>
#include <math.h>

/**
 *@brief XYZ - structure containing direction cosines coordinates (x, y, z)
 * x points East direction
 * y points North direction
 * z points Up
 **/

    typedef struct {
        float x, y, z;
    }equ2enu;


    typedef struct {
        float ha, dec;
        }hadec;


    equ2enu *XYZ_hadec(float x, float y, float z){
        equ2enu *pr = (equ2enu*)malloc(sizeof(equ2enu));
        pr->x = x;
        pr->y = y;
        pr->z = z;
        return pr;
    
        free(pr);
    }


/**
* @func that converts ha to enu direction cosines
 *
 *@In parameters
 *----------------
 *  ha - hour angle in raians.
 *  dec - declination angle  in radians.
 *  latitude - observers latitude in radians
 
 *@Out Parameters
 *----------------
 * structure containing direction cosines coordinates (x, y, z).
 **/


equ2enu *hadec_2_enu(hadec equ, float lat){
    float Xh,Yh, Zh;
    float ha = equ.ha;
    float dec = equ.dec;
    Xh = -sinf(ha) * cosf(dec);
    Yh = sinf(dec) * cosf(lat) - sinf(lat) * cosf(dec) * cosf(ha);
    Zh = cosf(lat) * cosf(dec) * cosf(ha) + sinf(lat) * sinf(dec);
    
    return XYZ_hadec(Xh, Yh, Zh);
    
}
    


    //hirax_enu *hadec_2_enu(hirax_hadec equ, float lat);

#endif
