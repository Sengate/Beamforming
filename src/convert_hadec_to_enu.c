
#include "convert_hadec_to_enu.h"

/*hirax_enu *hadec_2_enu(float ha, float dec, float lat){
    
    float Xh,Yh, Zh;
    
    Xh = -sinf(ha) * cosf(dec);
    Yh = sinf(dec) * cosf(lat) - sinf(lat) * cosf(dec) * cosf(ha);
    Zh = cosf(lat) * cosf(dec) * cosf(ha) + sinf(lat) * sinf(dec);
    
    return XYZ(Xh, Yh, Zh);
    
}*/


hirax_enu *hadec_2_enu(hirax_hadec equ, float lat){
    float Xh,Yh, Zh;
    float ha = equ.ha;
    float dec = equ.dec;
    Xh = -sinf(ha) * cosf(dec);
    Yh = sinf(dec) * cosf(lat) - sinf(lat) * cosf(dec) * cosf(ha);
    Zh = cosf(lat) * cosf(dec) * cosf(ha) + sinf(lat) * sinf(dec);
    
    return XYZ(Xh, Yh, Zh);
    
}
