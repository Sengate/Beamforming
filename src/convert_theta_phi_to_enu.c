#include "convert_theta_phi_to_enu.h"
#include <math.h>


//single precision
equ2enu *hadec_2_enu(hadec equ, float lat){
    float Xh,Yh, Zh;
    float ha = equ.ha;
    float dec = equ.dec;
    Xh = -sinf(ha) * cosf(dec);
    Yh = sinf(dec) * cosf(lat) - sinf(lat) * cosf(dec) * cosf(ha);
    Zh = cosf(lat) * cosf(dec) * cosf(ha) + sinf(lat) * sinf(dec);
    
    return XYZ_hadec(Xh, Yh, Zh);
    
}

