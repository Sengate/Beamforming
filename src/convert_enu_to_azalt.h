#ifndef __hirax_beamform_altaz
#define __hirax_beamform_altaz

#include <math.h>
    /**
     *@brief Spherical coorinates: azimuth, altitude in radians.
     **/
    typedef struct {
        float az;  /** The azimuth angle in radians **/
        float alt; /** The altitude angle in radians.**/
    }altaz;


    /**
     *@brief 3D Horizontal coordinates whose components are in .... units.
     **/
    typedef struct {
        float x;     /**The cartesian x-coordinate in .... units **/
        float y;     /**The cartesian x-coordinate in ....  units**/
        float z;     /**The cartesian x-coordinate in ....  units**/
    }enu_cords;

/*
 
    enu_cords *ENU_cords(float x, float y, float z){
        enu_cords *pr = (enu_cords*)malloc(sizeof(enu_cords));
        pr->x = x;
        pr->y = y;
        pr->z = z;
        return pr;
    
        free(pr);
    }
*/

    /**
     *@brief An function that convert Horizontal coordinates to spherical cordinates (az, alt).
     
     *@param[in]: struct enu containing;
        xh - the  east-west direction horizontal cordinates.
        yh - the  north-south direction horizontal cordinates.
        zh - the - up direction horizontal cordinates.
     
     *@param[out]: struct spherical coordinates (az, alt).
        az - azimuth angels in radians.
        alt - altitude angels in radians.
     **/

    altaz *hor2Altaz(equ2enu *u){
        
        float az = atan(u->x/u->y);
        float alt = acos(u->z);
        altaz *pe = (altaz*)malloc(sizeof(altaz));
        pe->az = az;
        pe->alt = alt;
        return pe;
        
        free(pe);
    }


#endif
