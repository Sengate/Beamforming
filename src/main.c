#include "convert_hadec_to_enu.h"
#include "convert_theta_phi_to_enu.h"
#include "convert_enu_to_azalt.h"
//#include "steering_vector.h"
#include "beams_steeringVec.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <omp.h>
#include <math.h>


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
            
            float angx = omega * (DirecCos->x * dat->xPositions[iant]);
            float angy = omega * (DirecCos->y * dat->yPositions[iant]);
            
            //east-west phases
            EWphases[ifreq][iant].re = cos(angx);
            EWphases[ifreq][iant].im = -1*sin(angx);
            
            //north-south phases
            NSphases[ifreq][iant].re = cos(angy);
            NSphases[ifreq][iant].im = -1*sin(angy);
            
            }
        }
    }




int main(int argc, char *argv[]){
    int i,j;
    float theta[] = {1.0,2.0,3.0};
    float phi[] = {1.0,2.0,5.0};
    int nbeam =3;
    float lat= 2.0;
    
    float xpos[] = {1,2,3};
    float ypos[] = {1,2,3};
    
  
    
    float freq[] ={0,1,2,3,4};
    
    //float (*out1)[5] = malloc(sizeof(out1[3][5]));
    //float out[3][5];
    
    Data_parameters *dat = (Data_parameters*)malloc(sizeof(Data_parameters));
    memset(dat, 0 ,sizeof(Data_parameters));
    
    dat->nchan = 5;
    dat->xPositions = xpos;
    dat->yPositions = ypos;
    dat->frequencies = freq;
    dat->nants =3;
    
    

   /* phase_memalloc(dat);
   
    steering_vector(dat, 1,2,EWphases, NSphases);
    
    for (unsigned int ifreq =0; ifreq < dat->nchan; ifreq++){
        for (unsigned int iant = 0; iant < dat->nants; iant++){
            
         // printf("%f ", NSphases[ifreq][iant].im);
            
        }
        
    }*/
    /*
    complex (*EW_beams_stv)[dat->nchan][dat->nants];
    
    make3Dmatrix_complex(nbeam, dat->nchan, dat->nants, &EW_beams_stv);
    
    
    for (i=0; i<3; i++){
        steering_vector(dat, theta[i], phi[i], EWphases, NSphases);
        for (unsigned int ifreq =0; ifreq < dat->nchan; ifreq++){
            for (unsigned int iant = 0; iant < dat->nants; iant++){
                EW_beams_stv[i][ifreq][iant].im = EWphases[ifreq][iant].re;
                
                printf("%f ", EW_beams_stv[i][ifreq][iant].im);
        //EW_beams_stv[i].re=EWphases.re;
            }}}*/
    //phase_destroy(dat);

    //steering_vector(dat, 1,2);
    //for (int i = 0; i<3; i++)
      //printf("%f ", dat->yPositions[i]*5);
    //for (int i = 0; i<3; i++)
     //for (int j = 0; j<5; j++)
     //  printf("%f ",out1[i][j]);
    
    
    
    free(dat);
    //phase_destroy(dat);
    //hirax_hadec *equ = (hirax_hadec*)malloc(sizeof(hirax_hadec));
    
   /* float *rr = malloc(sizeof(float) * 3);
    
    for (i=0; i<3; i++){
        hadec EDU = {theta[i], phi[i]};
        equ2enu *eq = hadec_2_enu(EDU, lat);
        altaz *out = hor2Altaz(eq);
        rr[i] = out->az;
        
         printf("%f ", out->az);
    }
    
    for (i=0; i<3; i++){
        printf("%f ", rr[i]);    }*/
    
    //equ2enu *eq = hadec_2_enu(EDU, lat);
        
        //altaz *out = hor2Altaz(eq);
        
        //printf("%f, ", out->alt);
        //equ2enu *rr = (equ2enu*)malloc(sizeof(equ2enu));
        
        //rr = hadec_2_enu(EDU, lat);
        //altaz *pp = hor2Altaz(rr);
        
        
        //thtph HOR ={ha[i], dec[i]};
        //printf("\n");
        //printf("%f", rr->x/rr->y);
        //float e = EDU.ha/HOR.theta;
        //printf("%f", e);
        //printf("\n");
        //xh[i] = EDU.ha;
        //yh[i]= EDU.ha;
        //zh[i] = EDU.dec * 3;
        
       //hirax_enu *rr = thtph_to_enu(EDU);
        //xh[i] = HADEC_2_enu->x;
        
        //printf("%f", HADEC_2_enu->x*100);
    
    //for (i=0; i<3; i++){
      //  printf("%f", zh[i]);
        
    //}
    
    //free(xh);
    //free(yh);
    //free(zh);
}
