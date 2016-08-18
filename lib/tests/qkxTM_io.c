#include <lime.h>

int read_custom_binary_gauge_field (double **gauge, char *fname, LatticeInfo *lattInfo)
{
/*
read gauge fileld config stored in binary file
*/
	FILE			*fid;
	int			x1, x2, x3, x4, ln[4] = { 0, 0, 0, 0 };
	unsigned long long	lvol, ixh, iy, mu;
	char			tmpVar[20];
	double			*U = (double*)malloc(18*sizeof(double));
	double			*resEvn[4], *resOdd[4];
	int			nvh, iDummy;
	LimeReader		*limereader;
	char			*lime_type, *lime_data;
	n_uint64_t		lime_data_size;
	char			dummy;
	int			isDouble;
	int			error_occured=0;
	double			dDummy;
	double			*ftmp=NULL;
	unsigned long long	iread, idx;


	fid=fopen(fname,"r");
	if(fid==NULL)
	{
		fprintf(stderr,"Error reading configuration! Could not open %s for reading\n",fname);
		error_occured=1;
	}
	
	if ((limereader = limeCreateReader(fid))==NULL)
	{
		fprintf(stderr,"Could not create limeReader\n");
		error_occured=1;
	}
	
	if(!error_occured)
	{
		while(limeReaderNextRecord(limereader) != LIME_EOF )
		{
			lime_type	= limeReaderType(limereader);

			if(strcmp(lime_type,"ildg-binary-data")==0)
				break;

			if(strcmp(lime_type,"xlf-info")==0)
			{
				lime_data_size = limeReaderBytes(limereader);
				lime_data = (char * )malloc(lime_data_size);
				limeReaderReadData((void *)lime_data,&lime_data_size, limereader);

				strcpy	(tmpVar, "kappa =");
				sscanf(qcd_getParamComma(tmpVar,lime_data, lime_data_size),"%lf",&dDummy);    
				printfQuda("Kappa:    \t%lf\n", dDummy);
				inv_param->kappa	= dDummy;

				strcpy	(tmpVar, "mu =");
				sscanf(qcd_getParamComma(tmpVar,lime_data, lime_data_size),"%lf",&dDummy);    
				printfQuda("Mu:       \t%lf\n", dDummy);

				if      (overrideMu)
				{
					printfQuda	("MU OVERRIDEN\t%lf\n", newMu);
					inv_param->mu	 = newMu;
				}
				else
					inv_param->mu		= dDummy;

				free(lime_data);
			}

			if(strcmp(lime_type,"ildg-format")==0)
			{
				lime_data_size = limeReaderBytes(limereader);
				lime_data = (char * )malloc(lime_data_size);
				limeReaderReadData((void *)lime_data,&lime_data_size, limereader);

				strcpy	(tmpVar, "<precision>");
				sscanf(qcd_getParam(tmpVar,lime_data, lime_data_size),"%i",&isDouble);    
				printfQuda("Precision:\t%i bit\n",isDouble);

				strcpy	(tmpVar, "<lx>");
				sscanf(qcd_getParam(tmpVar,lime_data, lime_data_size),"%i",&iDummy);
				param->X[0]	 = iDummy/gridSize[0];
				ln[0]		 = iDummy;

				strcpy	(tmpVar, "<ly>");
				sscanf(qcd_getParam(tmpVar,lime_data, lime_data_size),"%i",&iDummy);
				param->X[1]	 = iDummy/gridSize[1];
				ln[1]		 = iDummy;

				strcpy	(tmpVar, "<lz>");
				sscanf(qcd_getParam(tmpVar,lime_data, lime_data_size),"%i",&iDummy);
				param->X[2]	 = iDummy/gridSize[2];
				ln[2]		 = iDummy;

				strcpy	(tmpVar, "<lt>");
				sscanf(qcd_getParam(tmpVar,lime_data, lime_data_size),"%i",&iDummy);
				param->X[3]	 = iDummy/gridSize[3];
				ln[3]		 = iDummy;

				printfQuda("Volume:   \t%ix%ix%ix%i\n", ln[0], ln[1], ln[2], ln[3]);
				printfQuda("Subvolume:\t%ix%ix%ix%i\n", param->X[0], param->X[1], param->X[2], param->X[3]);

				free(lime_data);
			}
		}

		// Read 1 byte to set file-pointer to start of binary data


		lime_data_size=1;
		limeReaderReadData(&dummy,&lime_data_size,limereader);
		offset = ftell(fid)-1;

		limeDestroyReader(limereader);
	}     



	if(error_occured)
		return	1;

	if	(isDouble == 32)
	{
		printf	("Error: Unsupported precision %d bits\n", isDouble);
		return	1;
	}

	nvh = (param->X[0] * param->X[1] * param->X[2] * param->X[3]) / 2;

	for(int dir = 0; dir < 4; dir++)
	{
		resEvn[dir] = gauge[dir];
		resOdd[dir] = gauge[dir] + nvh*18;
	} 

	lvol = ln[0]*ln[1]*ln[2]*ln[3];

	if(lvol==0)
	{
		fprintf(stderr, "[] Error, zero volume\n");
		return(5);
	}
 

	strcpy	(tmpVar, "native");


	//load time-slice by time-slice:

	chunksize	 = 4*3*3*sizeof(qcd_complex_16);
	ftmp		 = (char*) malloc(((unsigned int) chunksize*nvh*2));

	if	(ftmp == NULL)
	{
		fprintf(stderr,"Error in qcd_getGaugeLime! Out of memory, couldn't alloc %u bytes\n", (unsigned int) (chunksize*nvh*2));
		return	1;
	}
	else
		printf	("%d bytes reserved for gauge fields (%d x %u)\n", chunksize*nvh*2, chunksize, nvh*2);

	if	(MPI_File_read_all(mpifid, ftmp, 4*3*3*2*nvh*2, MPI_DOUBLE, &status) == 1)
		printf  ("Error in MPI_File_read_all\n");
	
	if      (4*3*3*2*nvh*2*sizeof(double) > 2147483648)
	{
		printf  ("File too large. At least %lu processes are needed to read thils file properly.\n", (4*3*3*2*nvh*2*sizeof(double)/2147483648)+1);
		printf  ("If some results are wrong, try increasing the number of MPI processes.\n");
	}

	if	(!qcd_isBigEndian())
		qcd_swap_8	((double*) ftmp,2*4*3*3*nvh*2);
#else
	ftmp	 = (double*)malloc(lvol*72*sizeof(double));

	if	(ftmp == NULL)
	{
		fprintf	(stderr, "Error, could not alloc ftmp\n");
		return	6;
	}
 
	iread	 = fread	(ftmp, sizeof(double), 72*lvol, fid);

	if	(iread != 72*lvol)
	{
		fprintf(stderr, "Error, could not read proper amount of data\n");
		return	7;
	}

	fclose(fid);

	if	(!qcd_isBigEndian())      
        	qcd_swap_8	((double*) ftmp,72*lvol);
#endif

	// reconstruct gauge field
	// - assume index formula idx = (((t*LX+x)*LY+y)*LZ+z)*(4*3*2*2) + mu*(3*2*2) + 2*(3*u+c)+r
	//   with mu=0,1,2,3; u=0,1; c=0,1,2, r=0,1

	iy	 = 0;
	ixh	 = 0;

#ifdef	MPI_READCONF
	i	 = 0;
#endif

	int lx1 = 0, lx2 = 0, lx3 = 0, lx4 = 0;

	for(x4 = x4start; x4 < x4end; x4++) 
	{  // t
		for(x3 = x3start; x3 < x3end; x3++) 
		{  // z
			for(x2 = x2start; x2 < x2end; x2++) 
			{  // y
				for(x1 = x1start; x1 < x1end; x1++) 
				{  // x
					int oddBit	 = (x1+x2+x3+x4) & 1;

//					iy		 = x1+x2*ln[0]+x3*ln[0]*ln[1]+x4*ln[0]*ln[1]*ln[2];
#ifdef	MPI_READCONF
					iy		 = ((x1-x1start)+(x2-x2start)*param->X[0]+(x3-x3start)*param->X[1]*param->X[0]+(x4-x4start)*param->X[0]*param->X[1]*param->X[2])/2;
#else
					iy		 = x1+x2*param->X[0]+x3*param->X[1]*param->X[0]+x4*param->X[0]*param->X[1]*param->X[2];
#endif
					for(mu = 0; mu < 4; mu++)
					{  
#ifdef	MPI_READCONF
						memcpy(U, &(ftmp[i]), 18*sizeof(double));

						if(oddBit)
							memcpy(&(resOdd[mu][18*iy]), U, 18*sizeof(double));
						else
							memcpy(&(resEvn[mu][18*iy]), U, 18*sizeof(double));

						i	+= 144;
#else
						ixh		 = (lx1+lx2*param->X[0]+lx3*param->X[0]*param->X[1]+lx4*param->X[0]*param->X[1]*param->X[2])/2;


						idx = mu*18 + 72*iy;

						double *links_ptr = ftmp+idx;

						memcpy(U, links_ptr, 18*sizeof(double));

						if(oddBit)
							memcpy(resOdd[mu]+18*ixh, U, 18*sizeof(double));
						else
							memcpy(resEvn[mu]+18*ixh, U, 18*sizeof(double));
#endif
					}

					++lx1;
				}

				lx1	 = 0;
				++lx2;
			}

			lx2	 = 0;
			++lx3;				
		}

		lx3	 = 0;
		++lx4;
	}

	free(ftmp);

#ifdef	MPI_READCONF
	MPI_File_close(&mpifid);
#endif  

	//	Apply BC here

	applyGaugeFieldScaling<double>((double**)gauge, nvh, param);
  
	return	0;
}
