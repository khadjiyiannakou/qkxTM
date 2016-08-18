#include <stdio.h>
#include <stdlib.h>

int main(){

  FILE *ptr;
  ptr = fopen("kale.bin","rb");
  char *buffer = (char*)malloc(5*5*5*10*4*3*3*2*sizeof(double));
  int count = fread(buffer, sizeof(double), 5*5*5*10*4*3*3*2,ptr);
  printf("%d\n",count);
}
