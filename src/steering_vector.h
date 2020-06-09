#include "convert_hadec_to_enu.h"
#include "convert_theta_phi_to_enu.h"
#include "convert_enu_to_azalt.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <omp.h>
#include <math.h>

//define spped of light (c) and 2*pi
#define c 30000000.0
#define TWO_PI 3.14159265359


/**
 *@brief A structure of complex number.
 **/
typedef float real;
typedef struct {
    real re;
    real im;
}complex;


/**
 *@brief A structure containing antenna positions and frequencies.
 **/
typedef struct{
    int nants;           /* Number of antennas in 1 direction. */
    int nchan;           /* number of frequency channels. */
    float *xPositions;   /* Antenna positions in the east-west direction in ...units */
    float *yPositions;   /* Antenna positions in the nort-south direction in ... units */
    float *zPositions;   /* Antenna positions in the up direction  in .....*/
    float *frequencies;  /* Frequencies in .... */
}Data_parameters;


/* Complex antenna phases */
complex **EWphases; //east-west phases    dim: [nchan x nants].
complex **NSphases; //north-south phases  dim: [nchan x nants].



/**
 *@brief phase memory allocation.
 **/
void phase_memalloc(Data_parameters *dat){
    
    /* malloc east west phases */
    EWphases= (complex **)malloc(dat->nchan * sizeof(complex *));
    for (int i = 0; i<dat->nchan; i++)
        EWphases[i] = (complex *)malloc(dat->nants * sizeof(complex));
    
    /*malloc NS phases*/
    NSphases= (complex **)malloc(dat->nchan * sizeof(complex *));
    for (int i = 0; i<dat->nchan; i++)
        NSphases[i] = (complex *)malloc(dat->nants * sizeof(complex));
}

/**
 *@brief phase memory free.
 **/
void phase_destroy(Data_parameters *dat){
    
    for (int i = 0; i < dat->nchan; i++)
        free(EWphases[i]);
    
    free( EWphases);
    
    
    for (int i = 0; i < dat->nchan; i++)
        free( NSphases[i]);
    
    free(NSphases);
}


/**
 *@brief Function that computes antenna phases for each direction, and single frequency.
 
 *@param[in]:
 Data_parameters-structure contains Antenna positions and frequencies.
 theta - azimuth angle in radians.
 phi -  elevation angle in radians.
 *@param[out]:
 EWphases - complex phases for each frequency and antenna in the EW direction.
 NSphase - complex phases for each frequency and antenna in the NS direction.
 **/
void steering_vector(Data_parameters *dat, float theta, float phi, complex **EWphases, complex **NSphases){
    
    /*--------Euler's formula e^(i*ang) = cos(ang)-i*sin(ang)----------*/
    
    thtph tp = {theta, phi};
    thtph2enu *DirecCos = thtph_to_enu(tp);
    
    for (unsigned int ifreq =0; ifreq < dat->nchan; ifreq++){
        
        float omega =  TWO_PI * dat->frequencies[ifreq];
        
        for (unsigned int iant = 0; iant < dat->nants; iant++){
            
            float angx = omega * (DirecCos->x * dat->xPositions[iant])/c;
            float angy = omega * (DirecCos->y * dat->yPositions[iant])/c;
            
            //east-west phases
            EWphases[ifreq][iant].re = cos(angx);
            EWphases[ifreq][iant].im = -1*sin(angx);
            
            //north-south phases
            NSphases[ifreq][iant].re = cos(angy);
            NSphases[ifreq][iant].im = -1*sin(angy);
            
            //printf("%f ", EWphases[ifreq][iant].im);
            
        }
    }
}
