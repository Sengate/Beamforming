
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
