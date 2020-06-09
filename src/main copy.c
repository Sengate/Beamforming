


void make3Dmatrix_f(unsigned dim1, unsigned dim2, unsigned dim3, float (**ptr)[dim2][dim3]){
  *ptr = malloc(dim1 *sizeof(**ptr));
    
}


void  make4Dmatrix_f(unsigned dim1, unsigned dim2, unsigned dim3,  unsigned dim4, float (**ptr)[dim2][dim3][dim4]){
    *ptr = malloc(dim1 *sizeof(**ptr));
}


void  make5Dmatrix_f(unsigned dim1, unsigned dim2, unsigned dim3,  unsigned dim4,unsigned dim5, float (**ptr)[dim2][dim3][dim4][dim5]){
    *ptr = malloc(dim1 *sizeof(**ptr));
}
