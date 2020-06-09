#include "convert_enu_to_azalt.h"
#include <math.h>


//Get the azimuth and altitude angles from Horizontal coordinates (Xh, Yh, Zh)
altaz *hor2Altaz(equ2enu *u){
    float az = atan(u->x/u->y);
    float alt = acos(u->z);
    altaz *pe = (altaz*)malloc(sizeof(altaz));
    pe->az = az;
    pe->alt = alt;
    return pe;
    
    free(pe);
}
