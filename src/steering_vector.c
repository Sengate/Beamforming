#include "convert_hadec_to_enu.h"
#include "convert_theta_phi_to_enu.h"
#include "convert_enu_to_azalt.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <omp.h>
#include <math.h>

#define c 30000000.0
#define two_PI 2 * 3.14159265359


typedef float real;
typedef struct {
    real re;
    real im;
}complex


typedef struct{
    int nants;
    int nchan;
    float *xPositions;
    float *yPositions;
    float *zPositions;
    float *channel;
}Data;


typedef struct out_phases{
    float ** EWPhases;
    float ** NSPhases;
} out_phases;


void steering_vector(Data *dat, float theta, float phi){
    
    /*--------Euler's formula e^(i*ang) = cos(ang)+i*sin(ang)----------*/
    
    float *omega = (float*)malloc(sizeof(float)*dat->nchan);
    
    thtph tp = {theta, phi};
    thtph2enu *DirecCos = htph_to_enu(tp);
    
    for (unsigned int ifreq =0; ifreq < dat->nchan; ifreq++){
        
        omega[ifreq] = TWO_PI * dat->channel[ifreq];
        
        for (unsigned int iant = 0; iant < dat->iants; iant++){
            
            float angx = omega * (DirecCos->x * dat-xPositions[iant])/c;
            
            float angy = omega * (DirecCos->y * dat-yPositions[iant])/c;
            
            
            
            
    
    
    
}
