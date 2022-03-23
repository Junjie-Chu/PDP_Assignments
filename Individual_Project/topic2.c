#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>
#include <string.h>
// the prop provided by the teacher is used via prop.h
#include "prop.h"

// define the constant values
#define num_bins 20
#define total_time 100

// print the hists while size of data is small
void horizontal(int range[],int array[],int len){
    printf("\t↑\n");
    printf("\t|\n");
    printf("\t|\n");
    printf("\t|%d\n",range[0]);
    for(int i=0;i<len;i++){
       int k=0;   
       printf("\t|");
       for(;k<array[i];k++){
          printf("*");
       }
       printf("(%d)",array[i]);
       printf("\n\t|%d\n",range[i+1]);
    }
    printf("\t--------------------------------------------------->\n");
}

//main start
int main(int argc, char** argv){

	//check the number of input
	if(argc != 3){
		printf("Error:2 argvs needed! 1st argv is the size of data in each process; 2nd is the output file name \n"); 
		return -1; 
        }
        
        //initial the mpi communication world
	int size, rank;
    	MPI_Status status; 
    	MPI_Init(&argc, &argv); 
    	MPI_Comm_size(MPI_COMM_WORLD, &size); 
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    	
    	//variables we need in the simulation
    	int localsize = atoi(argv[1]);
    	int globalsize = localsize*size;
    	double T = total_time;
    	double t = 0;//current time
    	double a0,u1,u2,tau;
    	int r,max,min;
    	
    	//set seed of rand
    	srand(time(NULL)+rank);
    	//srand(rank);
    	
    	//static arrays we need in the simulation
    	int X[7],binrange[num_bins+1];
    	int X0[7] = {900, 900, 30, 330, 50, 270, 20};
    	double w[15];
    	//dynamic arrays we need in the simulation
	int *P;
	int *finalstate;
	int *samples;
	P = (int*)malloc((7*15)*sizeof(int));
	memset(P,0,(7*15)*sizeof(int));
	finalstate = (int*)malloc((7*localsize)*sizeof(int));
	memset(finalstate,0,(7*localsize)*sizeof(int));
	samples = (int*)malloc(localsize*sizeof(int));
	memset(samples,0,localsize*sizeof(int));
	//set the initial value of P
    	P[0*7+0] = 1; P[1*7+0] = -1; P[2*7+0] = -1; 
    	P[3*7+1] = 1; P[4*7+1] = -1; P[5*7+1] = -1;     
    	P[2*7+2] = 1; P[5*7+3] = 1; P[6*7+2] = -1; 
    	P[7*7+2] = -1; P[7*7+4] = 1; P[8*7+3] = -1; 
    	P[9*7+3] = -1; P[9*7+5] = 1; P[10*7+4] = -1; 
    	P[11*7+4] = -1; P[11*7+6] = 1; P[12*7+5] = -1; 
    	P[13*7+0] = 1; P[13*7+6] = -1; P[14*7+6] = -1;
    	
    	MPI_Barrier(MPI_COMM_WORLD);
    	double start = MPI_Wtime();
    	//start the simulation 	
	for(int i=0; i<localsize; i++){
		//step1:set T,t,X=X0
		//initial X<-X0
		for (int j = 0;j < 7;j++){
			X[j]=X0[j];
		}
		
		//start time <- 0
		t = 0;
		
		//step2:while condition:t<T
		while(t<T){
			//step3:w=prob(x)
			prop(X,w);
			
			//step4:a0 = sum w(i)-w(R)
	    		a0 = 0; 
	    		for(int j=0; j<15; j++){
				a0 += w[j]; 
	    		}
	    		//printf("\na0:%lf\n",a0);
	    		//step5:generate random numbers u1 and u2
	    		
	    		u1 = (double)rand()/RAND_MAX; 
	    		u2 = (double)rand()/RAND_MAX;
	    		//different accuracy
	    		//u1 = (float)rand()/RAND_MAX; 
	    		//u2 = (float)rand()/RAND_MAX;
	    		
	    		//step6:set tau = -ln(u1)/a0
	    		tau = -log(u1)/a0;
	    		
	    		//step7:find r
	    		double sumr = 0;
	    		for(int j=0; j<15; j++){
				sumr += w[j]; 
				if(sumr >= a0*u2){
	   				r = j;
	   				break;  
				}
    			} 
	    		//printf("\nr:%d\n",r);
	    		//step8:update X
	    		for(int j=0; j<7; j++){
				X[j] += P[r*7+j];
				//printf("\nX:%d\n",X[j]); 
	    		}
	    		
	    		//step9:update t
	    		t = t + tau; 
					
		}
		
		//save the final state
		memcpy(&(finalstate[i*7]), X, 7*sizeof(int)); 
	}
	
	//collect the susceptible samples in each process, the first element in each row of final state
	for (int j=0; j < localsize;j++){
		samples[j]=finalstate[j*7];
		//printf("\nsamples:%d\n",samples[j]);
	}
	
	//find the max and min of samples in each process
	//max
	max = 0;
	for (int j=0; j <localsize;j++){
		if(samples[j]>max){
			max = samples[j];
		}
	}
	//min
	min = samples[0];
	for (int j=0; j <localsize;j++){
		if(samples[j]<min){
			min = samples[j];
		}
	}
	//printf("\nmax:%d,min:%d\n",max,min);
	//summary the max and min in all processes
	MPI_Allreduce(&max, &(binrange[num_bins+1]), 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); 
    	MPI_Allreduce(&min, &(binrange[0]), 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD); 
	
	//find the global intervals(width of each bin),Upper bound(global max) and lower bound(global min)
	double width,globalmax,globalmin;
	globalmax = (double)binrange[num_bins+1];
	globalmin = (double)binrange[0];
	width = (double)((globalmax-globalmin)/num_bins);
	
	//set upper bound and lower bound of each interval
	for(int j=0; j<num_bins+1; j++){
		binrange[j] = (int)((int)globalmin+(int)(j*width));
		//printf("\nbinrange[%d]:%d\n",j,binrange[j]); 
    	}
    	
    	//according to the new range, collect the value of each bin
	int* binvalue;
	binvalue = (int*)malloc(num_bins*sizeof(int));
	memset(binvalue,0,num_bins*sizeof(int));
	int used = 0;
	for (int k=0;k<localsize;k++){
		int mark = 0;
		for(int j=0;j<num_bins;j++){
			if(binrange[j] <= samples[k] && samples[k] < binrange[j+1]){
				mark = 1;
				used++;
				binvalue[j]++;
				break;	
			}	
		}
		/*printf("\nused num:%d mark:%d\n",used,mark);
		if(mark==0){
			printf("I am samples[%d], my value is %d",k,samples[k]);
		}*/
	}
	//The max value havent be counted (max = binrange[20])
	if (localsize-used>0){
		//printf("localsize-used:%d",localsize-used);
		binvalue[num_bins-1] = binvalue[num_bins-1] + localsize-used;
	}
	
	//rank 0 collect the data
	int* globalbinvalue;
    	if(rank == 0){
    	/*
    		printf("\n");
    		for (int k=0;k<localsize;k++){
    			printf("%d,",samples[k]);
    		}
    		printf("\n");
    	*/
		globalbinvalue =(int*)malloc(sizeof(int)*num_bins);
		memset(globalbinvalue,0,sizeof(int)*num_bins);		 
    	}
    	MPI_Reduce(binvalue, globalbinvalue, num_bins, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);	
	
	//End of main process
	MPI_Barrier(MPI_COMM_WORLD);
	double timecost = MPI_Wtime() - start; 
	double global_time; 
	MPI_Reduce(&timecost, &global_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	
	//Print hists when size is small and save global bin values to file
	if(rank == 0){
		//time cost
		printf("Time cost: %f. \n", global_time);
		//print hist
		if (globalsize<=400){
			horizontal(binrange,globalbinvalue,num_bins);
		}
		
		//write output file
		FILE* output;
		int wstatus;
		output = fopen(argv[2],"w");
		for(int j=0; j<num_bins+1; j++){
			wstatus=fprintf(output, "%d, ", binrange[j]);
    		}
		fprintf(output, "\n");
		for(int j=0; j<num_bins; j++){
			wstatus=fprintf(output, "%d, ", globalbinvalue[j]);
    		}
		fprintf(output, "\n");
		fclose(output);
		free(globalbinvalue); 
    	}
    	
    	//free
    	free(P);
    	free(finalstate);
    	free(samples);
    	free(binvalue);
	
	//MPI END
	MPI_Finalize();
}
