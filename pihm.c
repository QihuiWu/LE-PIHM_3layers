 /* File        : pihm.c                                                        *
 * Version     : Nov, 2007 (2.0)                                               *
 * Developer of PIHM2.0:        Mukesh Kumar (muk139@psu.edu)                  *
 * Developer of PIHM1.0:        Yizhong Qu   (quyizhong@gmail.com)             *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 *..............MODIFICATIONS/ADDITIONS in PIHM 2.0............................*
 * a) All modifications in physical process representations  in this version   *
 *    are listed as header in f.c and is_sm_et.c.     			       *
 * b) All addition/modifications in variable and structure definition/declarat-*
 *    -ion are listed as header in read_alloc.c and initialize.c	       *
 * c) 3 new input files have been added for geology, landcover and calibration *
 *    data								       *
 * d) Ported to Sundials 2.1.0                                                 *
 *-----------------------------------------------------------------------------*
 *                                                                             *
 * PIHM is an integrated finite volume based hydrologic model. It simulates    *
 * channel routing, overland flow, groundwater flow, macropore based infiltra- *
 * tion and stormflow, throughfall, evaporation from overlandflow-subsurface-  *
 * canopy, transpiration and  snowmelt by full coupling of processes.          *
 * It uses semi-discrete finite volume approach to discretize PDEs into ODEs,  *
 * and henceforth solving the global system of ODEs using CVODE. Global ODEs   *
 * are created in f.c. Any modifications in the process equations has to be    *
 * performed in f.c
 *                                                                             *
 *-----------------------------------------------------------------------------*
 * For questions or comments, please contact                                   *
 *      --> Mukesh Kumar (muk139@psu.edu)                                      *
 *      --> Prof. Chris Duffy (cxd11@psu.edu)                                  *
 * This code is free for research purpose only.                                *
 * Please provide relevant references if you use this code in your research work*
 *-----------------------------------------------------------------------------*
 *									       *
 * DEVELOPMENT RELATED REFERENCES:					       *
 * PIHM2.0:								       *
 *	a) Kumar, M., 2008, "Development and Implementation of a Multiscale,   *
 *	Multiprocess Hydrologic Model". PhD Thesis, Penn State University      *
 *	b) Kumar, M, G.Bhatt & C.Duffy, "Coupling of Data and Processes in     *
 *	Mesoscale Watershed", Advances in Water Resources (submitted)          *
 * PIHM1.0:								       *
 *	a) Qu, Y., 2005, "An Integrated hydrologic model for multiproces       *
 *	simulation using semi-discrete finite volume approach".PhD Thesis, PSU *
 *	b) Qu, Y. & C. Duffy, 2007, "A semidiscrete finite volume formulation  *
 *	for multiprocess watershed simulation". Water Resources Research       *
 *-----------------------------------------------------------------------------*
 * LICENSE:
 *******************************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/* SUNDIAL Header Files */
#include "sundials_types.h"	/* realtype, integertype, booleantype
				 * defination */
#include "cvode.h"		/* CVODE header file                             */
#include "cvode_spgmr.h"	/* CVSPGMR linear header file                    */
#include "sundials_smalldense.h"/* use generic DENSE linear solver for
				 * "small"   */
#include "nvector_serial.h"	/* contains the definition of type N_Vector      */
#include "sundials_math.h"	/* contains UnitRoundoff, RSqrt, SQR
				 * functions   */
#include "cvode_dense.h"	/* CVDENSE header file                           */
#include "sundials_dense.h"	/* generic dense solver header file              */
#include "pihm.h"		/* Data Model and Variable Declarations     */
#define UNIT_C 1440		/* Unit Conversions */



/* Main Function */ 
int
main(int argc, char *argv[])
{
	Model_Data      mData;	/* Model Data                */
	Control_Data    cData;	/* Solver Control Data       */
	N_Vector        CV_Y;	/* State Variables Vector    */
	void           *cvode_mem,*bpdata;	/* Model Data Pointer        */
	N_Vector        CV_Y_Hydro;   /* State Variables Vector    */
    void           *cvode_mem_Hydro;      /* Model Data Pointer        */
	N_Vector        CV_Y_LE;   /* State Variables Vector    */
    void           *cvode_mem_LE;      /* Model Data Pointer        */
	int             flag;	/* flag to test return value */
	FILE           *Ofile[18];	/* Output file     */
	char           *ofn[18];
	FILE 		   *init1,*init2, *para_file,*newelev1,*newelev2;
	char           *init_char1,*init_char2,*para_file_char, *directory, *path,*temp_path,*newelev_char1,*newelev_char2;
	char 		   *init_file;
	int				IsExist=0;
	int 			OutputID=1;
	realtype		total_t,total_t_LE, P;
	#ifdef LE_PIHM
		FILE           *Ofile_LE[10];	/* Output file     */
		char           *ofn_LE[10];	
	#endif
	
	#ifdef LE_PIHM_SED
		
	#endif
	FILE           *iproj;	/* Project File */
	int             N;	/* Problem size              */
	int				N_Hydro,N_LE; /* Problem size of LE processes*/
	int             i, j, k, w;/* loop index                */
	realtype        NextPtr, StepSize, total_NextPtr,total_LE_NextPtr;	/* stress period & step size */
	clock_t         start, finish;	/* system clock at points    */
	char           *filename;
	realtype	T;
	int	q;
	time_t current_time;

	/* Project Input Name */
	if (argc != 2) {
		iproj = fopen("projectName.txt", "r");
		if (iproj == NULL) {
			printf("\t\nUsage ./pihm project_name");
			printf("\t\n         OR              ");
			printf("\t\nUsage ./pihm, and have a file in the current directory named projectName.txt with the project name in it");
			exit(0);
		} else {
			filename = (char *) malloc(15 * sizeof(char));
			fscanf(iproj, "%s", filename);
			fclose(iproj);
		}
	}
	else {
		/* get user specified file name in command line */
		filename = (char *) malloc(strlen(argv[1]) * sizeof(char));
		strcpy(filename, argv[1]);
	}
	/* Open Output Files */
	directory = (char *) malloc(1024 * sizeof(char));
	path = (char *) malloc(1024 * sizeof(char));
	temp_path = (char *) malloc(1024 * sizeof(char));
	
	strcpy(path, "output/");
	i = mkdir(path,0755);
	sprintf(temp_path,"%s%d/",path,OutputID);
	i = mkdir(temp_path,0755);
	while(i==-1)
	{
		free(temp_path);
		temp_path = (char *) malloc(1024 * sizeof(char));
		OutputID += 1;
		sprintf(temp_path,"%s%d/",path,OutputID);
		i = mkdir(temp_path,0755);
	}
	
	strcpy(directory,temp_path);
	
	copy_inputs(directory,filename);
	
	
	ofn[0] = (char *) malloc((strlen(filename) + 1024) * sizeof(char));
	strcpy(ofn[0], directory);
	strcat(ofn[0], filename);
	Ofile[0] = fopen(strcat(ofn[0], ".GW"), "w");
	ofn[1] = (char *) malloc((strlen(filename) + 1024) * sizeof(char));
	strcpy(ofn[1], directory);
	strcat(ofn[1], filename);
	Ofile[1] = fopen(strcat(ofn[1], ".surf"), "w");
	ofn[2] = (char *) malloc((strlen(filename) + 1024) * sizeof(char));
	strcpy(ofn[2], directory);
	strcat(ofn[2], filename);
	Ofile[2] = fopen(strcat(ofn[2], ".et0"), "w");
	ofn[3] = (char *) malloc((strlen(filename) + 1024) * sizeof(char));
	strcpy(ofn[3], directory);
	strcat(ofn[3], filename);
	Ofile[3] = fopen(strcat(ofn[3], ".et1"), "w");
	ofn[4] = (char *) malloc((strlen(filename) + 1024) * sizeof(char));
	strcpy(ofn[4], directory);
	strcat(ofn[4], filename);
	Ofile[4] = fopen(strcat(ofn[4], ".et2"), "w");
	ofn[5] = (char *) malloc((strlen(filename) + 1024) * sizeof(char));
	strcpy(ofn[5], directory);
	strcat(ofn[5], filename);
	Ofile[5] = fopen(strcat(ofn[5], ".is"), "w");
	ofn[6] = (char *) malloc((strlen(filename) + 1024) * sizeof(char));
	strcpy(ofn[6], directory);
	strcat(ofn[6], filename);
	Ofile[6] = fopen(strcat(ofn[6], ".snow"), "w");
	ofn[7] = (char *) malloc((strlen(filename) + 1024) * sizeof(char));
	strcpy(ofn[7], directory);
	strcat(ofn[7], filename);
	Ofile[7] = fopen(strcat(ofn[7], ".unsat"), "w");
	ofn[8] = (char *) malloc((strlen(filename) + 1024) * sizeof(char));
	strcpy(ofn[8], directory);
	strcat(ofn[8], filename);
	Ofile[8] = fopen(strcat(ofn[8], ".Rech"), "w");
	ofn[9] = (char *) malloc((strlen(filename) + 1024) * sizeof(char));
	strcpy(ofn[9], directory);
	strcat(ofn[9], filename);
	Ofile[9] = fopen(strcat(ofn[9], ".infil"), "w");
	ofn[10] = (char *) malloc((strlen(filename) + 1024) * sizeof(char));
	strcpy(ofn[10], directory);
	strcat(ofn[10], filename);
	Ofile[10]=fopen(strcat(ofn[10], ".prep"),"w");
	ofn[11] = (char *)malloc((strlen(filename)+1024)*sizeof(char));
	strcpy(ofn[11], directory);
	strcat(ofn[11], filename);
	Ofile[11]=fopen(strcat(ofn[11], ".Flux1"),"w");
	ofn[12] = (char *)malloc((strlen(filename)+1024)*sizeof(char));
	strcpy(ofn[12], directory);
	strcat(ofn[12], filename);
	Ofile[12]=fopen(strcat(ofn[12], ".Flux2"),"w");
	ofn[13] = (char *)malloc((strlen(filename)+1024)*sizeof(char));
	strcpy(ofn[13], directory);
	strcat(ofn[13], filename);
	Ofile[13]=fopen(strcat(ofn[13], ".Flux3"),"w");
	ofn[14] = (char *)malloc((strlen(filename)+1024)*sizeof(char));
	strcpy(ofn[14], directory);
	strcat(ofn[14], filename);
	Ofile[14]=fopen(strcat(ofn[14], ".SubFlux1"),"w");
	ofn[15] = (char *)malloc((strlen(filename)+1024)*sizeof(char));
	strcpy(ofn[15], directory);
	strcat(ofn[15], filename);
	Ofile[15]=fopen(strcat(ofn[15], ".SubFlux2"),"w");
	ofn[16] = (char *)malloc((strlen(filename)+1024)*sizeof(char));
	strcpy(ofn[16], directory);
	strcat(ofn[16], filename);
	Ofile[16]=fopen(strcat(ofn[16], ".SubFlux3"),"w");
	ofn[17] = (char *)malloc((strlen(filename)+1024)*sizeof(char));
	strcpy(ofn[17], directory);
	strcat(ofn[17], filename);
	Ofile[17]=fopen(strcat(ofn[17], ".NetPrep"),"w");

	
	init_char1 = (char *)malloc((strlen(filename)+1024)*sizeof(char));
	init_char2 = (char *)malloc((strlen(filename)+1024)*sizeof(char));
	strcpy(init_char1, directory);
	strcat(init_char1, filename);
	strcat(init_char1, ".init");
	strcpy(init_char2, filename);
	strcat(init_char2, ".init");
	para_file_char = (char *)malloc((strlen(filename)+1024)*sizeof(char));
	strcpy(para_file_char, filename);
	para_file_char = strcat(para_file_char, ".para");
	
	#ifdef LE_PIHM
		ofn_LE[0] = (char *)malloc((strlen(filename)+1024)*sizeof(char));
		strcpy(ofn_LE[0], directory);
		strcat(ofn_LE[0], filename);
		Ofile_LE[0]=fopen(strcat(ofn_LE[0], ".grdelev"),"w");
		ofn_LE[1] = (char *)malloc((strlen(filename)+1024)*sizeof(char));
		strcpy(ofn_LE[1], directory);
		strcat(ofn_LE[1], filename);
		Ofile_LE[1]=fopen(strcat(ofn_LE[1], ".bedelev"),"w");
		ofn_LE[2] = (char *)malloc((strlen(filename)+1024)*sizeof(char));
		strcpy(ofn_LE[2], directory);
		strcat(ofn_LE[2], filename);
		Ofile_LE[2]=fopen(strcat(ofn_LE[2], ".SoilSubFlux1"),"w");
		ofn_LE[3] = (char *)malloc((strlen(filename)+1024)*sizeof(char));
		strcpy(ofn_LE[3], directory);
		strcat(ofn_LE[3], filename);
		Ofile_LE[3]=fopen(strcat(ofn_LE[3], ".SoilSubFlux2"),"w");
		ofn_LE[4] = (char *)malloc((strlen(filename)+1024)*sizeof(char));
		strcpy(ofn_LE[4], directory);
		strcat(ofn_LE[4], filename);
		Ofile_LE[4]=fopen(strcat(ofn_LE[4], ".SoilSubFlux3"),"w");
		ofn_LE[5] = (char *)malloc((strlen(filename)+1024)*sizeof(char));
		strcpy(ofn_LE[5], directory);
		strcat(ofn_LE[5], filename);
		Ofile_LE[5]=fopen(strcat(ofn_LE[5], ".Weathering"),"w");
		ofn_LE[6] = (char *)malloc((strlen(filename)+1024)*sizeof(char));
		strcpy(ofn_LE[6], directory);
		strcat(ofn_LE[6], filename);
		Ofile_LE[6]=fopen(strcat(ofn_LE[6], ".Uplift"),"w");
		ofn_LE[7] = (char *)malloc((strlen(filename)+1024)*sizeof(char));
		strcpy(ofn_LE[7], directory);
		strcat(ofn_LE[7], filename);
		Ofile_LE[7]=fopen(strcat(ofn_LE[7], ".BedloadFlux1"),"w");
		ofn_LE[8] = (char *)malloc((strlen(filename)+1024)*sizeof(char));
		strcpy(ofn_LE[8], directory);
		strcat(ofn_LE[8], filename);
		Ofile_LE[8]=fopen(strcat(ofn_LE[8], ".BedloadFlux2"),"w");
		ofn_LE[9] = (char *)malloc((strlen(filename)+1024)*sizeof(char));
		strcpy(ofn_LE[9], directory);
		strcat(ofn_LE[9], filename);
		Ofile_LE[9]=fopen(strcat(ofn_LE[9], ".BedloadFlux3"),"w");
		
		newelev_char1 = (char *)malloc((strlen(filename)+1024)*sizeof(char));
		newelev_char2 = (char *)malloc((strlen(filename)+1024)*sizeof(char));
		strcpy(newelev_char1, directory);
		strcat(newelev_char1, filename);
		strcat(newelev_char1, ".newelev");
		strcpy(newelev_char2, filename);
		strcat(newelev_char2, ".newelev");
	#endif
		
	#ifdef LE_PIHM_SED
			
	#endif
	
	
	
	
	
	/* allocate memory for model data structure */
	mData = (Model_Data) malloc(sizeof *mData);

	printf("\n ...  PIHM 2.2 is starting ... \n");

	/* read in 9 input files with "filename" as prefix */
	read_alloc(filename, mData, &cData);

	/*
	 * if(mData->UnsatMode ==1) {    }
	 */
	#ifdef FULLY_COUPLE
		if (mData->UnsatMode == 2) 
		{
			/* problem size */
			N = 3 * mData->NumEle;
			#ifndef LE_PIHM_SUBHYDRO
				N = mData->NumEle;
			#endif
			#ifdef LE_PIHM
				N = 5 * mData->NumEle;
				#ifndef LE_PIHM_SUBHYDRO
					N = 3 * mData->NumEle;
				#endif
			#endif
			mData->DummyY = (realtype *) malloc(N * sizeof(realtype));
			/* initial state variable depending on machine */
			CV_Y = N_VNew_Serial(N);
			/* initialize mode data structure */
			initialize(filename, mData, &cData, CV_Y);
			/* allocate memory for solver */
			cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
			if (cvode_mem == NULL) 
			{
				printf("CVodeMalloc failed. \n");
				return (1);
			}
			flag = CVodeSetFdata(cvode_mem, mData);
			flag = CVodeSetInitStep(cvode_mem, cData.InitStep);
			flag = CVodeSetStabLimDet(cvode_mem, TRUE);
			flag = CVodeSetMaxStep(cvode_mem, cData.MaxStep);
			flag = CVodeMalloc(cvode_mem, f, cData.StartTime, CV_Y, CV_SS, cData.reltol, &cData.abstol);
			flag = CVSpgmr(cvode_mem, PREC_NONE, 0);
			
			// flag = CVodeSetFdata(cvode_mem, mData);
			// flag = CVodeMalloc(cvode_mem, f, cData.StartTime, CV_Y, CV_SS, cData.reltol, &cData.abstol);
			// bpdata = CVBandPrecAlloc (cvode_mem, N, 3, 3);
			// flag = CVBPSpgmr(cvode_mem, PREC_LEFT, 0, bpdata);
		}
	#endif
	#ifndef FULLY_COUPLE
		if (mData->UnsatMode == 2) 
		{
			/* problem size */
			N_Hydro = 3 * mData->NumEle;
			#ifndef LE_PIHM_SUBHYDRO
				N_Hydro = mData->NumEle;
			#endif
			N_LE = 2 * mData->NumEle;
			mData->DummyY_Hydro = (realtype *) malloc(N_Hydro * sizeof(realtype));
			mData->DummyY_LE = (realtype *) malloc(N_LE * sizeof(realtype));
			/* initial state variable depending on machine */
			CV_Y_Hydro = N_VNew_Serial(N_Hydro);
			CV_Y_LE = N_VNew_Serial(N_LE);
			/* initialize mode data structure */
			initialize_no_fully(filename, mData, &cData, CV_Y_Hydro, CV_Y_LE);
			/* allocate memory for solver */
			cvode_mem_Hydro = CVodeCreate(CV_BDF, CV_NEWTON);
			if (cvode_mem_Hydro == NULL) 
			{
				printf("CVodeMalloc failed. \n");
				return (1);
			}
			flag = CVodeSetFdata(cvode_mem_Hydro, mData);
			flag = CVodeSetInitStep(cvode_mem_Hydro, cData.InitStep);
			flag = CVodeSetStabLimDet(cvode_mem_Hydro, TRUE);
			flag = CVodeSetMaxStep(cvode_mem_Hydro, cData.MaxStep);
			flag = CVodeMalloc(cvode_mem_Hydro, f_Hydro, cData.StartTime, CV_Y_Hydro, CV_SS, cData.reltol, &cData.abstol);
			flag = CVSpgmr(cvode_mem_Hydro, PREC_NONE, 0);
			P = NV_Ith_S(CV_Y_Hydro, 0);
			
			cvode_mem_LE = CVodeCreate(CV_BDF, CV_NEWTON);
			if (cvode_mem_LE == NULL) 
			{
				printf("CVodeMalloc failed. \n");
				return (1);
			}
			flag = CVodeSetFdata(cvode_mem_LE, mData);
			flag = CVodeSetInitStep(cvode_mem_LE, cData.InitStep);
			flag = CVodeSetStabLimDet(cvode_mem_LE, TRUE);
			flag = CVodeSetMaxStep(cvode_mem_LE, cData.MaxStep);
			flag = CVodeMalloc(cvode_mem_LE, f_LE, cData.StartTime, CV_Y_LE, CV_SS, cData.reltol, &cData.abstol);
			flag = CVSpgmr(cvode_mem_LE, PREC_NONE, 0);
		}
	#endif	
	FILE *Filename[3];
	Filename[0] = fopen("coordinate_x.txt","w");
	Filename[1] = fopen("coordinate_y.txt","w");
	Filename[2] = fopen("elevation.txt","w");
  	for ( k=0;k<mData->NumEle;k++)
	{
		fprintf(Filename[0],"%lf\t",mData->Ele[k].x);
               fprintf(Filename[0],"\n"); 
		fprintf(Filename[1],"%lf\t",mData->Ele[k].y);
               fprintf(Filename[1],"\n");
		fprintf(Filename[2],"%lf\t",mData->Ele[k].zmax);
               fprintf(Filename[2],"\n");         
	}
	fflush(Filename[0]); 
	fflush(Filename[1]); 
	fflush(Filename[2]); 
	fclose(Filename[0]);
	fclose(Filename[1]);
	fclose(Filename[2]);
	
	printf("\nSolving ODE system ... \n");
	
	// start = clock();

	/* start solver in loops */
	total_t = cData.StartTime;
	total_t_LE = cData.StartTime;
	while (total_t < cData.FinishTime)
	{
		/* set start time */
		mData->t = cData.StartTime;
		for (i = 0; i < cData.NumSteps; i++) 
		{
		/*
		 * inner loops to next output points with ET step size
		 * control
		 */
			while (mData->t < cData.Tout[i + 1]) {
				if (mData->t + cData.ETStep >= cData.Tout[i + 1]) {
					NextPtr = cData.Tout[i + 1];
				} else {
					NextPtr = mData->t + cData.ETStep;
				}
				StepSize = NextPtr - mData->t;
				total_NextPtr = total_t + StepSize;
				total_LE_NextPtr = total_t_LE + StepSize;
				/* calculate Interception Storage */
				
				printf("\n Tsteps = %f	total t = %f", mData->t, total_t);
				
				#ifdef FULLY_COUPLE
					if (mData->t>9)
					{
						printf("gogogo");
					}
					is_sm_et(mData->t, StepSize, mData, CV_Y);
					flag = CVode(cvode_mem, total_NextPtr, CV_Y, &total_t, CV_NORMAL);
				#else
					is_sm_et(mData->t, StepSize, mData, CV_Y_Hydro);
					flag = CVode(cvode_mem_Hydro, total_NextPtr, CV_Y_Hydro, &total_t, CV_NORMAL);
					for(k=0; k<mData->NumEle; k++)
					{
						if (NV_Ith_S(CV_Y_Hydro, k)<0) NV_Ith_S(CV_Y_Hydro, k)=0;
						if (NV_Ith_S(CV_Y_Hydro, k + mData->NumEle)<0) NV_Ith_S(CV_Y_Hydro, k + mData->NumEle)=0;
						if (NV_Ith_S(CV_Y_Hydro, k + 2*mData->NumEle)<0) NV_Ith_S(CV_Y_Hydro, k + 2*mData->NumEle)=0;
						
						mData->SurfDepth[k] = NV_Ith_S(CV_Y_Hydro, k); 
						mData->GW[k] = NV_Ith_S(CV_Y_Hydro, k + mData->NumEle);
						mData->UGW[k] = NV_Ith_S(CV_Y_Hydro, k + 2*mData->NumEle);
					} 
					flag = CVode(cvode_mem_LE, total_LE_NextPtr, CV_Y_LE, &total_t_LE, CV_NORMAL);
					for(k=0; k<mData->NumEle; k++)
					{
						mData->grdelev[k] = NV_Ith_S(CV_Y_LE, k); 
						mData->bedelev[k] = NV_Ith_S(CV_Y_LE, k + mData->NumEle);
					}
				#endif
				
				mData->t = NextPtr;			
				update(mData->t, mData);
			}
		
			#ifdef FULLY_COUPLE
				PrintData(Ofile, &cData, mData, CV_Y, total_t);
				#ifdef LE_PIHM
					PrintData_LE(Ofile_LE, &cData, mData, CV_Y, total_t);
				#endif				
			#else
				PrintData_no_fully(Ofile, &cData, mData, CV_Y_Hydro, total_t);
				PrintData_no_fully_LE(Ofile_LE, &cData, mData, CV_Y_LE, total_t_LE);
			#endif
			
			if ((int) total_t%cData.groundelevDInt == 0)
			{
				//printf("\n Tsteps = %f	total t = %f", mData->t, total_t);
				init1 = fopen(init_char1,"w");
				newelev1 = fopen(newelev_char1,"w");
				
				init2 = fopen(init_char2,"w");
				newelev2 = fopen(newelev_char2,"w");
				
				para_file = fopen(para_file_char,"w"); 
				
				for (k=0;k<mData->NumEle;k++)
				{
					#ifdef FULLY_COUPLE
						#ifdef LE_PIHM
							#ifdef LE_PIHM_SUBHYDRO
								fprintf(init1,"%lf\t%lf\t%lf\t%lf\t%lf",mData->EleIS[k], mData->EleSnowCanopy[k]+mData->EleSnowGrnd[k], NV_Ith_S(CV_Y, k),NV_Ith_S(CV_Y, k+mData->NumEle),NV_Ith_S(CV_Y, k+2*mData->NumEle));
									   fprintf(init1,"\n"); 
									   
								fprintf(newelev1,"%lf\t%lf",NV_Ith_S(CV_Y, k+3*mData->NumEle),NV_Ith_S(CV_Y, k+4*mData->NumEle));
									   fprintf(newelev1,"\n");
								
								fprintf(init2,"%lf\t%lf\t%lf\t%lf\t%lf",mData->EleIS[k], mData->EleSnowCanopy[k]+mData->EleSnowGrnd[k], NV_Ith_S(CV_Y, k),NV_Ith_S(CV_Y, k+mData->NumEle),NV_Ith_S(CV_Y, k+2*mData->NumEle));
									   fprintf(init2,"\n"); 
									   
								fprintf(newelev2,"%lf\t%lf",NV_Ith_S(CV_Y, k+3*mData->NumEle),NV_Ith_S(CV_Y, k+4*mData->NumEle));
									   fprintf(newelev2,"\n");
							#else
								fprintf(init1,"%lf\t%lf\t%lf\t0\t0",mData->EleIS[k], mData->EleSnowCanopy[k]+mData->EleSnowGrnd[k], NV_Ith_S(CV_Y, k));
									   fprintf(init1,"\n"); 
									   
								fprintf(newelev1,"%lf\t%lf",NV_Ith_S(CV_Y, k+mData->NumEle),NV_Ith_S(CV_Y, k+2*mData->NumEle));
									   fprintf(newelev1,"\n");
								
								fprintf(init2,"%lf\t%lf\t%lf\t0\t0",mData->EleIS[k], mData->EleSnowCanopy[k]+mData->EleSnowGrnd[k], NV_Ith_S(CV_Y, k));
									   fprintf(init2,"\n"); 
									   
								fprintf(newelev2,"%lf\t%lf",NV_Ith_S(CV_Y, k+mData->NumEle),NV_Ith_S(CV_Y, k+2*mData->NumEle));
									   fprintf(newelev2,"\n");
							#endif
						#else
							#ifdef LE_PIHM_SUBHYDRO
								fprintf(init1,"%lf\t%lf\t%lf\t%lf\t%lf",mData->EleIS[k], mData->EleSnowCanopy[k]+mData->EleSnowGrnd[k], NV_Ith_S(CV_Y, k),NV_Ith_S(CV_Y, k+mData->NumEle),NV_Ith_S(CV_Y, k+2*mData->NumEle));
									   fprintf(init1,"\n"); 
									   
								fprintf(newelev1,"%lf\t%lf",mData->Ele[k].zmax,mData->Ele[k].zmin);
									   fprintf(newelev1,"\n");
								
								fprintf(init2,"%lf\t%lf\t%lf\t%lf\t%lf",mData->EleIS[k], mData->EleSnowCanopy[k]+mData->EleSnowGrnd[k], NV_Ith_S(CV_Y, k),NV_Ith_S(CV_Y, k+mData->NumEle),NV_Ith_S(CV_Y, k+2*mData->NumEle));
									   fprintf(init2,"\n"); 
									   
								fprintf(newelev2,"%lf\t%lf",mData->Ele[k].zmax,mData->Ele[k].zmin);
									   fprintf(newelev2,"\n");
							#else
								fprintf(init1,"%lf\t%lf\t%lf\t0\t0",mData->EleIS[k], mData->EleSnowCanopy[k]+mData->EleSnowGrnd[k], NV_Ith_S(CV_Y, k));
									   fprintf(init1,"\n"); 
									   
								fprintf(newelev1,"%lf\t%lf",mData->Ele[k].zmax,mData->Ele[k].zmin);
									   fprintf(newelev1,"\n");
								
								fprintf(init2,"%lf\t%lf\t%lf\t0\t0",mData->EleIS[k], mData->EleSnowCanopy[k]+mData->EleSnowGrnd[k], NV_Ith_S(CV_Y, k));
									   fprintf(init2,"\n"); 
									   
								fprintf(newelev2,"%lf\t%lf",mData->Ele[k].zmax,mData->Ele[k].zmin);
									   fprintf(newelev2,"\n");
							#endif
						#endif
					#else
						fprintf(init1,"%lf\t%lf\t%lf\t%lf\t%lf",mData->EleIS[k], mData->EleSnowCanopy[k]+mData->EleSnowGrnd[k], NV_Ith_S(CV_Y_Hydro, k),NV_Ith_S(CV_Y_Hydro, k+mData->NumEle),NV_Ith_S(CV_Y_Hydro, k+2*mData->NumEle));
								fprintf(init1,"\n"); 
								   
						fprintf(newelev1,"%lf\t%lf",NV_Ith_S(CV_Y_LE, k),NV_Ith_S(CV_Y_LE, k+mData->NumEle));
								fprintf(newelev1,"\n");
							
						fprintf(init2,"%lf\t%lf\t%lf\t%lf\t%lf",mData->EleIS[k], mData->EleSnowCanopy[k]+mData->EleSnowGrnd[k], NV_Ith_S(CV_Y_Hydro, k),NV_Ith_S(CV_Y_Hydro, k+mData->NumEle),NV_Ith_S(CV_Y_Hydro, k+2*mData->NumEle));
								fprintf(init2,"\n"); 
								   
						fprintf(newelev2,"%lf\t%lf",NV_Ith_S(CV_Y_LE, k),NV_Ith_S(CV_Y_LE, k+mData->NumEle));
								fprintf(newelev2,"\n");
					#endif
				}
				
				fprintf(para_file, "Verbose\t%d\n", cData.Verbose);
                fprintf(para_file, "Debug\t%d\n",cData.Debug);
                fprintf(para_file, "Hydro_init_type\t2\n");
                fprintf(para_file, "LE_init_type\t1\n");
                fprintf(para_file, "groundelevDInt\t%d\n",cData.groundelevDInt);
                fprintf(para_file, "rockelevDInt\t%d\n",cData.rockelevDInt);
                fprintf(para_file, "bedloadDInt\t%d\n",cData.bedloadDInt);
                fprintf(para_file, "soilsubDInt\t%d\n",cData.soilsubDInt);
                fprintf(para_file, "upliftDInt\t%d\n",cData.upliftDInt);
                fprintf(para_file, "weatheringDInt\t%d\n",cData.weatheringDInt);
                fprintf(para_file, "gwDInt\t%d\n",cData.gwDInt);
                fprintf(para_file, "surfDInt\t%d\n",cData.surfDInt);    
                fprintf(para_file, "snowDInt\t%d\n",cData.snowDInt);
                fprintf(para_file, "prepDInt\t%d\n",cData.prepDInt);    
                fprintf(para_file, "RechInt\t%d\n",cData.RechInt);
                fprintf(para_file, "IsDInt\t%d\n",cData.IsDInt);
                fprintf(para_file, "usDInt\t%d\n",cData.usDInt);
                fprintf(para_file, "etInt\t%d\n",cData.etInt);
                fprintf(para_file, "computtimeDInt\t%d\n",cData.computtimeDInt);
                fprintf(para_file, "UnsatMode\t%d\n",mData->UnsatMode);
                fprintf(para_file, "SurfMode\t%d\n",mData->SurfMode);
                fprintf(para_file, "Solver\t%d\n",cData.Solver);
                fprintf(para_file, "GSType\t%d\n",cData.GSType);
                fprintf(para_file, "MaxK\t%d\n",cData.MaxK);
                fprintf(para_file, "delt\t%f\n",cData.delt);
                fprintf(para_file, "abstol\t%f\n",cData.abstol);
                fprintf(para_file, "reltol\t%f\n",cData.reltol);
                fprintf(para_file, "InitStep\t%f\n",cData.InitStep);
                fprintf(para_file, "MaxStep\t%f\n",cData.MaxStep);
                fprintf(para_file, "ETStep\t%f\n",cData.ETStep);
                if (mData->t<cData.EndTime) 
                {
                    fprintf(para_file, "StartTime\t%f\n",mData->t);
                }else
                {
                    fprintf(para_file, "StartTime\t%f\n",cData.StartTime);
                }
                fprintf(para_file, "EndTime\t%f\n",cData.EndTime);
				fprintf(para_file, "FinishTime\t%f\n",cData.FinishTime);
                fprintf(para_file, "outtype\t%d\n",cData.outtype);
                fprintf(para_file, "MSF\t%lf\n",mData->MSF);
                fprintf(para_file, "a\t%f\n",cData.a);
                fprintf(para_file, "b\t%f\n",cData.b);              
				
				
				
				fflush(init1); 
				fflush(newelev1); 
				fclose(init1);
				fclose(newelev1);
				fflush(init2); 
				fflush(newelev2); 
				fclose(init2);
				fclose(newelev2);
				fflush(para_file);
				fclose(para_file);
			}
		}
		for(k=0; k<mData->NumPrep; k++)
  		{
			mData->TSD_Prep[k].iCounter=0; 
  		}  
  
		for(k=0; k<mData->NumTemp; k++)
  		{   		
    		mData->TSD_Temp[k].iCounter=0;    		
  		} 
  
		for(k=0; k<mData->NumHumidity; k++)
  		{    		
    		mData->TSD_Humidity[k].iCounter=0;   		
  		} 
  
		for(k=0; k<mData->NumWindVel; k++)
  		{
    		mData->TSD_WindVel[k].iCounter=0;   		
  		} 

		for(k=0; k<mData->NumRn; k++)
  		{
    		mData->TSD_Rn[k].iCounter=0;   		
  		} 

		for(k=0; k<mData->NumG; k++)
  		{
    		mData->TSD_G[k].iCounter=0;    		
  		} 

		for(k=0; k<mData->NumP; k++)
  		{
    		mData->TSD_Pressure[k].iCounter=0;    		
  		} 

		for(k=0; k<mData->NumLC; k++)
  		{
    		mData->TSD_LAI[k].iCounter=0;   		
  		} 

		for(k=0; k<mData->NumLC; k++)
  		{
    		mData->TSD_RL[k].iCounter=0;
		}			
  
		for(k=0; k<mData->NumMeltF; k++)
  		{
    		mData->TSD_MeltF[k].iCounter=0;    		
  		}  		
	}

	
	/* Free memory */
	
	#ifndef FULLY_COUPLE
		N_VDestroy_Serial(CV_Y_Hydro);
		N_VDestroy_Serial(CV_Y_LE);
		CVodeFree(&cvode_mem_Hydro);
		CVodeFree(&cvode_mem_LE);
	#else
		N_VDestroy_Serial(CV_Y);
		CVodeFree(&cvode_mem);
	#endif
	//N_VDestroy_Serial(CV_Ydot);
	/* Free integrator memory */
	
	FreeData(mData, &cData);
	
	    for(i=0;i<18;i++)fclose(Ofile[i]);
		
		#ifdef LE_PIHM
			for(i=0;i<10;i++)fclose(Ofile_LE[i]);
		#endif
		
        for(i=0;i<18;i++)free(ofn[i]);
		
        #ifdef LE_PIHM
			for(i=0;i<10;i++)free(ofn_LE[i]);
		#endif
		
        free(filename);
		free(mData);
		free(init_char1);
		free(newelev_char1);
		free(init_char2);
		free(newelev_char2);
		free(directory);
        return 0;
}
