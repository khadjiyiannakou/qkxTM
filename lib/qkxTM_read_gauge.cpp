#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <lime.h>
#include <qkxTM.h>
#include <unistd.h>


/*
  Function that returns 1 if machine is big endian,
  0 if little endina. (by Giannis Koutsou)
*/
 
int qkxTM_isBigEndian()
{
   union{
     char C[4];
     int  R   ;
        }word;
   word.R=1;
   if(word.C[3]==1) return 1;
   if(word.C[0]==1) return 0;
   return -1;
}


/*
  Function that converts between
  big endian and little endian
  8 byte words.
*/
void qkxTM_swap_8(double *Rd, int N)
{
   register char *i,*j,*k;
   char swap;
   char *max;
   char *R = (char*) Rd;

   max = R+(N<<3);
   for(i=R;i<max;i+=8)
   {
      j=i; k=j+7;
      swap = *j; *j = *k;  *k = swap;
      j++; k--;
      swap = *j; *j = *k;  *k = swap;
      j++; k--;
      swap = *j; *j = *k;  *k = swap;
      j++; k--;
      swap = *j; *j = *k;  *k = swap;
   }
}


void qkxTM_swap_4(float *Rd, int N)
{
  register char *i,*j,*k;
  char swap;
  char *max;
  char *R =(char*) Rd;

  max = R+(N<<2);
  for(i=R;i<max;i+=4)
  {
    j=i; k=j+3;
    swap = *j; *j = *k;  *k = swap;
    j++; k--;
    swap = *j; *j = *k;  *k = swap;
  }
}


/* by Giannis Koutsou */
char* qkxTM_getParams(char* fname,int *len)
{
   FILE *pfile;
   char *params;
   int i;
   size_t size;

   if ((pfile=fopen(fname,"r"))==NULL)
   {
      fprintf(stderr,"Error, cannot open %s for reading\n",fname);
      return(NULL);
   }
  
   i=0;
   while(!feof(pfile))
   {
      fgetc(pfile);
      i++;
   }
   *(len)=i;
   rewind(pfile);
   params = (char*)malloc(*(len)*sizeof(char));

   size = fread(params,sizeof(char),*len,pfile);
   
   fclose(pfile);
  
   return params;
}


/* by Giannis Koutsou */
char* qkxTM_getParam(char token[],char* params,int len)
{
   int i,token_len=strlen(token);

   for(i=0;i<len-token_len;i++)
   {
      if(memcmp(token,params+i,token_len)==0)
      {
         i+=token_len;
         *(strchr(params+i,'<'))='\0';
         break;
      }
   }
   return params+i;
}
 
char* qkxTM_getParamComma(char token[],char* params,int len)
{
  int i,token_len=strlen(token);

  for(i=0;i<len-token_len;i++)
    {
      if(memcmp(token,params+i,token_len)==0)
	{
	  i+=token_len;
	  *(strchr(params+i,','))='\0';
	  break;
	}
    }
  return params+i;
}  



/* reads the text message from the header of ILDG configuration files */
int qkxTM_getLimeMessage(char *fname)
{
   LimeReader *limereader;
   FILE *fid;
   char *lime_type;//,*lime_data;
   n_uint64_t lime_data_size=0;
   int error_occured=0;
   char *message = NULL;
         
       /* read lime header */
   fid=fopen(fname,"r");
   if(fid==NULL)
     {
       fprintf(stderr,"Error in qkxTM_getLimeMessage! Could not open %s for reading\n",fname);
       error_occured=1;
       exit(EXIT_FAILURE);
     }else{
     if((limereader = limeCreateReader(fid))==NULL)
       {
         fprintf(stderr,"Error in qkxTM_getLimeMessage! Could not create limeReader\n");
         error_occured=1;
	 exit(EXIT_FAILURE);
       }
   }

   if(!error_occured)
     {
       while(limeReaderNextRecord(limereader) != LIME_EOF)
         {
	   lime_type = limeReaderType(limereader);
	   if(strcmp(lime_type,"ildg-binary-data")==0)
	     {
               break;
	     }
	   if(strcmp(lime_type,"xlf-info")==0)
	     {
               lime_data_size = limeReaderBytes(limereader);
               message = (char * )malloc(lime_data_size+1);
               limeReaderReadData((void *)message,&lime_data_size, limereader);
               message[lime_data_size]='\0';
               printf("lime-message:%s",message);
	     }
         }
       fclose(fid);
     }

   return(0);
}//end qkxTM_getLimeGmessage


 
int qkxTM_getGaugeLime(char *fname, qkxTMComplex **u, LatticeInfo *latInfo)
{
   LimeReader *limereader;
   FILE *fid;
   char *lime_type = NULL,*lime_data = NULL;
   n_uint64_t lime_data_size;
   char dummy;
   unsigned long long int offset;

   unsigned int chunksize;
   char *buffer;
   unsigned int  isDouble;
   int lx,ly,lz,lt;

   double mu;
   double kappa;
   double epsilon = 1e-8;

   //   int error_occured=0;
      
   if(u == NULL)
   {
      fprintf(stderr,"Error in qkxTM_getGaugeLime guage field pointer not initialized properly \n");
      exit(EXIT_FAILURE);
   }   

   
   // read lime header 
   fid=fopen(fname,"r");
   if(fid==NULL)
     {
       fprintf(stderr,"Error in qkxTM_getGaugeLime! Could not open %s for reading\n",fname);
       exit(EXIT_FAILURE);
     }

   if ( (limereader = limeCreateReader(fid)) == NULL )
     {
       fprintf(stderr,"Error in qkxTM_getGaugeLime! Could not create limeReader\n");
       exit(EXIT_FAILURE);
     }

   
   
   while(limeReaderNextRecord(limereader) != LIME_EOF )  // read all lime headers
     {
       
       lime_type = limeReaderType(limereader);
       if(strcmp(lime_type,"ildg-binary-data")==0)
	 {
	   break;
	 }

       if(strcmp(lime_type,"ildg-format")==0)         // this if will give info about lattice size and precision
	 {
	   

	   lime_data_size = limeReaderBytes(limereader);
	   lime_data = (char * )malloc(lime_data_size);
	   limeReaderReadData((void *)lime_data,&lime_data_size, limereader);
	   
	   sscanf(qkxTM_getParam((char *)"<precision>",lime_data, lime_data_size),"%i",&isDouble);    
	   if(isDouble != 64){
	     fprintf(stderr,"Error support only double precision\n");
	     exit(EXIT_FAILURE);
	   }

	   sscanf(qkxTM_getParam((char *)"<lx>",lime_data, lime_data_size),"%i",&lx);
	   if(lx != latInfo->L[0]){
	     fprintf(stderr,"Error lattice dimensions must agree\n");
             exit(EXIT_FAILURE);
           }

	   sscanf(qkxTM_getParam((char *)"<ly>",lime_data, lime_data_size),"%i",&ly);
	   if(ly != latInfo->L[1]){
	     fprintf(stderr,"Error lattice dimensions must agree\n");
             exit(EXIT_FAILURE);
           }

	   sscanf(qkxTM_getParam((char *)"<lz>",lime_data, lime_data_size),"%i",&lz);
	   if(lz != latInfo->L[2]){
	     fprintf(stderr,"Error lattice dimensions must agree\n");
             exit(EXIT_FAILURE);
           }

	   sscanf(qkxTM_getParam((char *)"<lt>",lime_data, lime_data_size),"%i",&lt);
	   if(lt != latInfo->L[3]){
	     fprintf(stderr,"Error lattice dimensions must agree\n");
             exit(EXIT_FAILURE);
           }

	   free(lime_data);
	 }

       if(strcmp(lime_type,"xlf-info")==0)
	 {
	   lime_data_size = limeReaderBytes(limereader);
	   lime_data = (char * )malloc(lime_data_size);
	   limeReaderReadData((void *)lime_data,&lime_data_size, limereader);


	   sscanf(qkxTM_getParamComma((char *)"kappa =",lime_data, lime_data_size),"%lf",&kappa);
	   if( (fabs(kappa - latInfo->kappa ) > epsilon) ){
	     fprintf(stderr," Error kappa values must agree\n");
             exit(EXIT_FAILURE);
           }
	  

	   sscanf(qkxTM_getParamComma((char *)"mu =",lime_data, lime_data_size),"%lf",&mu);
	     if( (fabs(mu - latInfo->mu ) > epsilon) ){
	     fprintf(stderr," Error mu values must agree\n");
             exit(EXIT_FAILURE);
           }
	   
	   
	   
	   free(lime_data);

	 }

     }



   // read 1 byte to set file-pointer to start of binary data 
   lime_data_size=1;
   limeReaderReadData(&dummy,&lime_data_size,limereader);   // i think this put the file pointer to the binary record
   offset = ftell(fid)-1;
   limeDestroyReader(limereader);      // i dont need limereader anymore
   fseek ( fid , offset , SEEK_SET );

   //   fclose(fid);
     
   
   //load time-slice by time-slice:


   chunksize=4*3*3*2*sizeof(double);
   unsigned long long int volume = lx * ly * lz * lt;
   buffer = (char*) malloc( chunksize * volume );
   size_t iread;

   if(buffer==NULL)
     {
       fprintf(stderr,"Error in qkxTM_getGaugeLime! Out of memory\n");
       exit(EXIT_FAILURE);
     }

   //   MPI_File_read_all(mpifid, buffer, 4*3*3*2*u->geo->lV, MPI_DOUBLE, &status);
   iread = fread(buffer, sizeof(double), 4*3*3*2*volume , fid);

   fprintf(stdout,"From conf %d number of double has read\n",(int)iread);

   if(!qkxTM_isBigEndian()){      
     qkxTM_swap_8((double*) buffer,2*4*3*3*volume);
     printf("big endian\n");
   }
   //   i=0;

   //      qkxTMComplex **buffer_c = (qkxTMComplex**) buffer;

   //   double test;

   //for(int i = 0 ; i < 1000 ; i+=8){
   //test = *( (double*) &(buffer[i]) );
   //printf("%e\n",test);
   // }
   //	int lx1 = 0, lx2 = 0, lx3 = 0, lx4 = 0;

   long long int i=0;
   
      for(int x4 = 0; x4 < lt; x4++) // t
	for(int x3 = 0; x3 < lz; x3++) // z
	  for(int x2 = 0; x2 < ly; x2++) // y
	    for(int x1 = 0; x1 < lx; x1++) { // x
	    
	      int iv = x4*lz*ly*lx + x3*ly*lx + x2*lx + x1;
			
	      for(int mu = 0; mu < 4; mu++)
		for(int ic1 = 0 ; ic1 < 3 ; ic1++)
		  for(int ic2 = 0 ; ic2 < 3 ; ic2++){
		    u[mu][MG(iv,ic1,ic2)] = *( (qkxTMComplex*) &(buffer[i]) );
		    i+= (int) sizeof(qkxTMComplex);
		  }
	      //		    memcpy(u[mu], buffer_c + iv*4*3*3 + mu*3*3 , 9*sizeof(qkxTMComplex));
	      
	    }
   

   free(buffer);
   //   MPI_File_close(&mpifid);
   // MPI_Type_free(&subblock);
   fclose(fid);
   return 0;
}//end qkxTM_getGaugeLime 


int qkxTM_getGaugeAndreas(char *fname, qkxTMComplex **u, LatticeInfo *latInfo)
{
   FILE *fid; 

   unsigned int chunksize;
   char *buffer;

   int lx = latInfo->L[0];
   int ly = latInfo->L[1];
   int lz = latInfo->L[2];
   int lt = latInfo->L[3];
   //   int error_occured=0;
      
   if(u == NULL)
     {
       fprintf(stderr,"Error in qkxTM_getGaugeLime guage field pointer not initialized properly \n");
       exit(EXIT_FAILURE);
     }   

   
   // read lime header 
   fid=fopen(fname,"rb");
   if(fid==NULL)
     {
       fprintf(stderr,"Error in qkxTM_getGaugeLime! Could not open %s for reading\n",fname);
       exit(EXIT_FAILURE);
     }
   

   //   int offset = ftell(fid)-1;
   // fseek ( fid , offset , SEEK_SET );

   chunksize=4*3*3*2*sizeof(double);
   unsigned long long int volume = lx *ly *lz *lt;
   buffer = (char*) malloc( chunksize * volume );
   size_t iread;

   if(buffer==NULL)
     {
       fprintf(stderr,"Error in qkxTM_getGauge! Out of memory\n");
       exit(EXIT_FAILURE);
     }

   //   MPI_File_read_all(mpifid, buffer, 4*3*3*2*u->geo->lV, MPI_DOUBLE, &status);
   iread = fread(buffer, sizeof(double), 4*3*3*2*volume , fid);

   //   fprintf(stdout,"From conf %d number of double has read\n",(int)iread);

   if(!qkxTM_isBigEndian()){      
     qkxTM_swap_8((double*) buffer,2*4*3*3*volume);
     //     printf("big endian\n");
   }

   long long int i=0;
   

   for(int ic1 = 0 ; ic1 < 3 ; ic1++)
     for(int ic2 = 0 ; ic2 < 3 ; ic2++)		 
       for(int x1 = 0; x1 < lx; x1++) 
	 for(int x2 = 0; x2 < ly; x2++) 
	   for(int x3 = 0; x3 < lz; x3++) 
	     for(int x4 = 0; x4 < lt; x4++) { 
	       
	       int iv = x4*lz*ly*lx + x3*ly*lx + x2*lx + x1;
	     
	       for(int mu = 0; mu < 4; mu++){
		 u[mu][MG(iv,ic1,ic2)] = *( (qkxTMComplex*) &(buffer[i]) );
		 i+= (int) sizeof(qkxTMComplex);
	       }
	     //		    memcpy(u[mu], buffer_c + iv*4*3*3 + mu*3*3 , 9*sizeof(qkxTMComplex));
	     
	   }
   
   
   free(buffer);
   //   MPI_File_close(&mpifid);
   // MPI_Type_free(&subblock);
   fclose(fid);
   return 0;
}//end qkxTM_getGaugeLime 
