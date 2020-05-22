#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include <libgen.h>


#define speed_light 300000000.0
#define two_PI 2 * 3.14159265359

typedef struct{
    /*Antenna Positions*/
    int nelements;
    float *xPositions;
    float *yPositions;
    float *zPositions;
    
    /*TABs directions*/
    int nbeams;
    float *beams_RAs;
    float *beams_DECs;
    
    /*LST list*/
    int ntimes;
    float *LST_list;
    
    /*Frequencies*/
    int nchan;
    float *frequencies; /*!! file/calculate?*/
    
}Ingredients;


/*FUNC 1!*/
/*Read files*/
Ingredients *read_files(char *antPositionsFile, char *beamsDircFile, char *lstFile, int nbeams, int nant, int nlst, int nchan){
    

    Ingredients  *INPUTS = (Ingredients*)malloc(sizeof(Ingredients));
    
    INPUTS->nelements = nant;
    INPUTS->nchan = nchan;
    INPUTS->ntimes = nlst;
    INPUTS->nbeams = nbeams;
    
    
    
    int x,y,z;
    int ra, dec, lst;
    
    FILE *antPtr, *beamPtr, *lstPtr;
    
    /*---------- Read antenna positions file------------------------*/
    
    if ((antPtr = fopen(antPositionsFile, "r")) == NULL)
    {
        printf("Error : Unable to open Antenna positions file\n");
        exit(EXIT_FAILURE);
    }
    
    INPUTS->xPositions = (float*)malloc(sizeof(float) * INPUTS->nelements);
    INPUTS->yPositions = (float*)malloc(sizeof(float) * INPUTS->nelements);
    INPUTS->zPositions = (float*)malloc(sizeof(float) * INPUTS->nelements);
    
    int count =0;
    while ( ( fscanf(antPtr, "%d %d %d", &x, &y, &z) ) == 3)
    {
        INPUTS->xPositions[count] = x;
        INPUTS->yPositions[count] = y;
        INPUTS->zPositions[count] = z;
        
        count++;
        if ( count >= INPUTS->nelements)
        {
            break;
        }
    }
    
    /*-------------Read beams directions file-----------*/
    
    if ((beamPtr = fopen(beamsDircFile, "r")) == NULL)
    {
        printf("Error : Unable to open Beams direction file\n");
        exit(EXIT_FAILURE);
    }
    INPUTS->beams_RAs = (float*)malloc(sizeof(float) * INPUTS->nbeams);
    INPUTS->beams_DECs = (float*)malloc(sizeof(float) * INPUTS->nbeams);
    int count2 =0;
    while ( ( fscanf(beamPtr, "%d %d", &ra, &dec) ) == 2)
    {
        INPUTS->beams_RAs[count2] = ra;
        INPUTS->beams_DECs[count2] = dec;
        count2++;
        if ( count2 >= INPUTS->nbeams)
        {
            break;
        }
    }
    
    
    /*-------------------Read LST file--------------------*/
   
    if ((lstPtr = fopen(lstFile, "r")) == NULL)
    {
        printf("Error : Unable to open LST file\n");
        exit(EXIT_FAILURE);
    }
    INPUTS->LST_list = (float*)malloc(sizeof(float) * INPUTS->ntimes);
    
    
    int count3 =0;
    while ( ( fscanf(lstPtr, "%d", &lst) ) == 1)
    {
        INPUTS->LST_list[count3] = lst;
   
        count3++;
        if ( count3 >= INPUTS->ntimes)
        {
            break;
        }
    }
    
   //Free allocated memories
    free(INPUTS->xPositions);
    free(INPUTS->yPositions);
    free(INPUTS->zPositions);
    free(INPUTS->beams_RAs);
    free(INPUTS->beams_DECs);
    free(INPUTS->LST_list);
    free(INPUTS);
    return INPUTS;
    }



/*---------------------------------Coordinates Coversion-------------------------------------------------*/
typedef struct {
    float az, alt;
}AzAlt;

typedef struct {
    float l, n, m;
}Dircos;


Dircos *XYZ(float x, float y, float z){
    Dircos *pr = (Dircos*)malloc(sizeof(Dircos));
    pr->l = x;
    pr->m = y;
    pr->n = z;
    return pr;
    
    free(pr);
}

Dircos *SphrCord(AzAlt pe){
    float azimuth = pe.az;
    float altitude = pe.alt;
    float X = cosf(azimuth) * sinf(altitude);
    float Y = sinf(azimuth) * sinf(altitude);
    float Z = cosf(altitude);
    return XYZ(X, Y, Z);
}



Dircos *Equtorial_2_Horizontal(float ha, float dec, float lat){
    float Xh= sinf(dec) * cosf(lat) - sinf(lat) * cosf(dec) * cosf(ha);
    float Yh = -sinf(ha) * cosf(dec);
    float Zh = cosf(lat) * cosf(dec) * cosf(ha) + sinf(lat) * sinf(dec);
    return XYZ(Xh, Yh, Zh);
}


float get_ha(float lst, float ra){
    float ha = lst -ra;
    return ha;
}



AzAlt *Get_azalt(Dircos *equ2har){
    //Ang azalt;
    float az = atan(equ2har->m/equ2har->l);
    float alt = asin(equ2har->n);
    AzAlt *pe = (AzAlt*)malloc(sizeof(AzAlt));
    pe->az = az;
    pe->alt = alt;
    return pe;
    
    free(pe);
}

/*------------------------------Calculate Steering vectors---------------------------------------------*/
typedef float real;
typedef struct {
    real re;
    real im;
}complex;


void mem_allocate(Ingredients *Inputs){
    /* allocate steering vectors memory (2D)*/
    complex (*EW_phases)[Inputs->nelements] = malloc( sizeof(EW_phases[Inputs->nchan][Inputs->nelements]));
    complex (*NS_phases)[Inputs->nelements] = malloc( sizeof(NS_phases[Inputs->nchan][Inputs->nelements]));
    //complex (*NS_phases)[Inputs->nelements] = malloc(Inputs->nchan * sizeof(*NS_phases));
}


/*void mem_destroy(){
free(omega);
    free(EW_phases);
    free(NS_phases);
}*/



/*FUNC2*/
/*zenith steering vector*/

void compute_zenith_steeringVectors(Ingredients *Inputs, complex **EW_phases, complex **NS_phases, float azimuth, float altitude){
    
    /*--------Euler's formula e^(i*ang) = cos(ang)+i*sin(ang)----------*/
    float *omega = (float*)malloc(sizeof(float)*Inputs->nchan);
    /*get direction cosines*/
    AzAlt azalt = {azimuth, altitude};
    Dircos *spherical_coordinates = SphrCord(azalt);
    
    for (unsigned int i =0; i<Inputs->nchan; i++){
        
        omega[i] = two_PI * Inputs->frequencies[i];
        
        /* steering vector */
        for (unsigned int k=0; k<Inputs->nelements; k++){
            
            float angx = omega[i] * (spherical_coordinates->l * Inputs->xPositions[k])/speed_light;
            float angy = omega[i] * (spherical_coordinates->m * Inputs->yPositions[k])/speed_light;
            
            /* East West steering vector*/
            EW_phases[i][k].re = cos(angx);
            EW_phases[i][k].im = -1 * sin(angx);
            
            /* North south steering vector*/
            NS_phases[i][k].re = cos(angy);
            NS_phases[i][k].im = -1 * sin(angy);
            
        }}}

/*--------------------- Compute beams steering vectors----------------------*/


void  make_4D_matrix(int dim1, int dim2, int dim3, int dim4, complex (**ptr)[dim2][dim3][dim4]){
    *ptr = malloc(dim1 *sizeof(**ptr));
}

void beam_mem_allocattor(Ingredients * Inputs){
    complex (*EW_beams_strVec)[Inputs->nbeams][Inputs->nchan][Inputs->nelements];
    complex (*NS_beams_strVec)[Inputs->nbeams][Inputs->nchan][Inputs->nelements];
    
    make_4D_matrix(Inputs->ntimes,Inputs->nbeams, Inputs->nchan, Inputs->nelements, &EW_beams_strVec);
    make_4D_matrix(Inputs->ntimes,Inputs->nbeams, Inputs->nchan, Inputs->nelements, &NS_beams_strVec);
    
    float *ha = (float*)malloc(sizeof(float)*Inputs->ntimes);
}

/*void beams_memo_destoy(){
free(EW_beams_strVec);
free(NS_beams_strVec);
    free(ha);}*/


/*--------------------------------------------------------------------------*/
/*FUNC3*/
/*beams phases*/
void compute_beams_phases(Ingredients *Inputs, complex ****EW_beamsStv, complex ****NS_beamsStv, float lat ){
    
    //float *ha = (float*)malloc(sizeof(float)*Inputs->ntimes);
    
    
    for (int i=0; i<Inputs->ntimes; i++){
        
        for (int j=0; j<Inputs->nbeams; j++{
            
            ha[i] = Inputs->LST_list[i] - Inputs->beams_RAs[j];
            float dec = Inputs->beams_DECs[j];
            
            /*convert equatorial coordiates to horizontal*/
            Dircos *equ2hor=Equtorial_2_Horizontal(ha, dec, lat);
            
            /* get beams azimuths and latitudes */
            AzAlt *get_azalt = Get_azalt(beam_sphcord);
            AzAlt beams_azalt = {az_b, alt_b};
            
            /*get lmn*/
            Dircos *beam_sphcord = SphrCord(beams_azalt);
            
            
            for (int k=0; k<Inputs->nchan; k++){
        
                for(int n=0; n<Inputs->nelements; n++){
                    
                    float omega = two_PI * Inputs->frequencies[k];
                    
                    float angx = omega * (beam_sphcord->l * Inputs->xPositions[n])/speed_light;
                    float angy = omega * (beam_sphcord->m * Inputs->yPositions[n])/speed_light;
                    
                    EW_beamsStv[i][j][k][n].re = cos(angx);
                    EW_beamsStv[i][j][k][n].im = -sin(angy);
                    
                    NS_beamsStv[i][j][k][n].re = cos(angx);
                    NS_beamsStv[i][j][k][n].im = -sin(angy);
                    
                }
            }
            
        }
    }
             
}

 
/*--------------------------------------*/
/*             TO DO  LIST              */
/*--------------------------------------*/
             
/* Write to files
 * Test
 *Tracking-frequency alignment, tracking space
 */



int main(){
    
    
}
    
    /*
    Coord *out = (Coord*)malloc(sizeof(Coord));
    
    int i;
    float X;
    int Y;
    int Z;
    char *filename = "file1.txt";
    FILE *ptr;
    
    if ((ptr = fopen(filename, "r")) == NULL)
    {
        printf("Error : Unable to open file for reading\n");
        exit(EXIT_FAILURE);
    }
    int count =0;
    out->x = (float *)malloc(sizeof(float) * 60);
    out->y = malloc(sizeof(float) * 60);
    out->z = malloc(sizeof(float) * 60);
    //int size =numl(fileName);
    // printf("%d", size);
    // printf("%d", line);
    
    while ( ( fscanf(ptr, "%f %d %d", &X, &Y, &Z) ) == 3)
    {
        out->x[count] = X;
        out->y[count] = Y;
        out->z[count] = Z;
        
        count++;
        if ( count >= 60)
        {
            break;
        }
    }
    
 
    
    for (i=0;i<60; i++){
        printf("\n%f",out->x[i] );
    }
    free(out->x);
    free(out->y);
    free(out->z);
    
    free(out);
    */
  
    

