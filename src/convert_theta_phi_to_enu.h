
#ifndef __hirax_beamform_enu
#define __hirax_beamform_enu

//includes
#include <stdlib.h>
#include <math.h>

/*
 * @struct thtph - structure containing angles theta, phi
 
 *theta zenith angle in raians
 *phi azimuth angle, measured from East(x) to North(y), in radians.

*@struct XYZ - structure containing direction cosines coordinates (x, y, z)
 * x points East direction
 * y points North direction
 * z points Up
 */

typedef struct {
    float theta, phi;
}thtph;


typedef struct {
    float x, y, z;
}thtph2enu;


thtph2enu *XYZ_thtph(float x, float y, float z){
    thtph2enu *pr = (thtph2enu*)malloc(sizeof(thtph2enu));
    pr->x = x;
    pr->y = y;
    pr->z = z;
    return pr;
    
    free(pr);
}


/*
 * @func that converts spherical coordinates to enu direction cosines
 *
 *@In parameters
 *----------------
 *  theta  in raians.
 *  phi  in radians.
 
 *@Out Parameters
 *----------------
 * structure containing direction cosines coordinates (x, y, z).
 */

//single precision
thtph2enu *thtph_to_enu(thtph theta_phi){
    
    float theta = theta_phi.theta;
    float phi = theta_phi.phi;
    
    float x = cosf(theta) * sinf(phi);
    float y = sinf(theta) * sinf(phi);
    float z = cosf(theta);
    
    return XYZ_thtph(x, y, z);
}


//hirax_enu *thtph_to_enu(thtph theta_phi);

#endif
