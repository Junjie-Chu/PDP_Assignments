#define PI 3.14159265358979323846
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv) {
	if (3 != argc) {
		printf("Usage: stencil num_values num_steps\n");
		return 0;
	}
    int num_values = atoi(argv[1]);
    int num_steps = atoi(argv[2]);
    int myid;
    int size;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);

    int send_count;
    send_count = num_values / size;
    double h = 2.0*PI/num_values;
    //give extra 2+2 elements to local input and output for computing
    double *localinput=(double*)malloc((send_count+4)*sizeof(double));
    double *localoutput=(double*)malloc((send_count+4)*sizeof(double));
    // Stencil values
    const int STENCIL_WIDTH = 5;
    const int EXTENT = STENCIL_WIDTH/2;
    const double STENCIL[] = {1.0/(12*h), -8.0/(12*h), 0.0, 8.0/(12*h), -1.0/(12*h)}; 
    double* input;
    double* output;
    double CommS,CommE,StenS,StenE,TotalS,TotalE;
    TotalS = MPI_Wtime();
    // 0 process: generate data and scatter
    if (myid == 0){
    	printf("\nnumber of processes: %d \n",size);
    	printf("number of data in each process: %d \n",send_count);
    	// Generate values for stencil operation
	input=(double*)malloc(num_values*sizeof(double));
	//printf("\ntotal input:\n");
    	for (int i=0; i<num_values; i++) {
    		input[i]=sin(h*i);
    		//printf("%f,",input[i]);
    	}
    	// Allocate data for result
    	output =(double*)malloc(num_values*sizeof(double));
    	}
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Scatter the Values
    MPI_Scatter(input, send_count, MPI_DOUBLE, localinput+2, send_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    /*
    //check each initial input
    if (myid == 0){
    	printf("\nmy id: %d, localinput:\n",myid);
    	for (int i=0; i<send_count;i++){
    		printf("%f ",localinput[2+i]);
    	}
    }*/
    
    
    // Define left and right process
    int left, right;
    left = myid -1;
    if (left < 0){left = size -1;}
    right = myid + 1;
    if (right > size -1){right = 0;}
    //printf("\nmy id:%d left: %d right: %d\n", myid, left, right);
    
    // Start timer
    double start = MPI_Wtime();
    
    // Repeatedly apply stencil
    for (int s=0; s<num_steps; s++) {	
	// The process has send_count+4 elements in total. 
	// EXTRA:The first 2 is localinput/localinput+1. The last 2 is localinput+sendcount+3/localinput+sendcount+2. 
	// The useful part: from localinput+2 to localinput+sendcount+1.
	//CommS = MPI_Wtime();
	MPI_Sendrecv(localinput+2, 2, MPI_DOUBLE, left, 1, localinput+send_count+2, 2, MPI_DOUBLE, right, 1, MPI_COMM_WORLD, &status);
	MPI_Sendrecv(localinput+send_count, 2, MPI_DOUBLE, right, 0, localinput, 2, MPI_DOUBLE, left, 0, MPI_COMM_WORLD, &status);
	//CommE = MPI_Wtime();
	/*printf("\nbefore sendrecv\n");
	printf("\nafter sendrecv\n");
       
       // Check the big array after sendrecv
    	if (myid == 0){
    		printf("\nmy id: %d, localinput+4:\n",myid);
    		for (int i=0; i<send_count+4;i++){
    			printf("%f ",localinput[i]);
    		}
    	}*/
    		
    	//StenS = MPI_Wtime();
    	// Apply stencil on points
	for (int i=2; i<send_count+2; i++) {
		double result = 0;
		for (int j=0; j<STENCIL_WIDTH; j++) {
			int index = i - EXTENT + j;
			result += STENCIL[j] * localinput[index];
		}
		localoutput[i] = result;
	}
	//StenE = MPI_Wtime();
	/*
	// Apply stencil on left boundary with periodic cond
	for (int i=2; i<EXTENT+2; i++) {
		double result = 0;
		for (int j=0; j<STENCIL_WIDTH; j++) {
			int index = (i - EXTENT + j + send_count) % send_count;
			result += STENCIL[j] * localinput[index];
		}
		localoutput[i] = result;
	}
    
    	// Apply stencil on right boundary with periodic cond
	for (int i=send_count-EXTENT+2; i<send_count+2; i++) {
		double result = 0;
		for (int j=0; j<STENCIL_WIDTH; j++) {
			int index = (i - EXTENT + j) % send_count;
			if(myid==0){printf("\nindex: %d\n",index);}
			result += STENCIL[j] * localinput[index];
		}
		localoutput[i] = result;
	}
	*/
	/*
	// Check the big localoutput
    	if (myid == 0){
    		printf("\nmy id: %d, localoutput+4:\n",myid);
    		for (int i=0; i<send_count+4;i++){
    			printf("%f ",localoutput[i]);
    		}
    	}*/
    	
	// Swap input and output
	if (s < num_steps-1) {
		double *tmp = localinput;
		localinput = localoutput;
		localoutput = tmp;
	}
    	MPI_Barrier(MPI_COMM_WORLD);
    }
	
    // Stop timer
    double my_execution_time = MPI_Wtime() - start;
    double max_execution_time;
    MPI_Reduce( &my_execution_time, &max_execution_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
       
    MPI_Barrier(MPI_COMM_WORLD);
    //Gather data
    MPI_Gather(localoutput+2,send_count,MPI_DOUBLE,output,send_count,MPI_DOUBLE,0,MPI_COMM_WORLD);
    TotalE = MPI_Wtime();
    
    if(myid==0){
    	//printf("\nCommunication time in last timestep: %f\n", CommE-CommS);
    	//printf("\nStencil time in last timestep: %f\n", StenE-StenS);
    	printf("\nTotal_stencil_max_execution_time: %f\n", max_execution_time);
    	printf("\nTotal time: %f\n", TotalE-TotalS);
    } 
    
    //0 process write the .txt
    if (myid == 0){
     	//Write to file
     	FILE *file=fopen("outputP.txt","w");
     	for (int i = 0; i < num_values; i++)
         	fprintf(file, "%.4f \n", output[i]);
     	fclose(file);
     	
     	free(input);
    	free(output);
    }
    
    // Clean up
    free(localinput);
    free(localoutput);
    
    // End MPI
    MPI_Finalize();
    return 0;
}

