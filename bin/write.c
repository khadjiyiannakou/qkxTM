#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(){

  double *buffer = (double*)malloc(2000000*sizeof(double));
  memset(buffer,0,2000000*sizeof(double));
  FILE *ptr;
  ptr = fopen("kale.bin","wb");
  fwrite(buffer,2000000,sizeof(double),ptr);
  
}
