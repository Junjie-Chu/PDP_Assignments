#define PI 3.14159265358979323846
#define _XOPEN_SOURCE
#define STRAT_ONE ((int)1)
#define STRAT_TWO ((int)2)
#define STRAT_THREE ((int)3)
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>

int pivot_strat;
MPI_Status status;
int currentsize;

/* Forward Declaration */
double* gen_num_vector(int seq, int len);
double* merge(double *v1, int n1, double *v2, int n2);
int cmpfunc_desc(const void * a, const void * b);
int cmpfunc (const void * a, const void * b);
double find_pivot(double* data, int len);
double find_mean(double* data, int len);
void check_pivot(double* data, int len, double pivot, int rank);
double pivot(int pivot_strat,int len, double *rcv_buffer, int rank, int size, MPI_Comm comm);
double* parallel_quicksort(MPI_Comm group, int rank, int size, double* data, int chunk);

int main(int argc, char **argv) {

    int size, rank,len, seq, chunk, remainder;
    double *data, *rcv_buffer, l_pivot, parr_runtime, loc_runtime;
    /* Initializing MPI */
    MPI_Status status;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &size);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    seq=atoi(argv[1]);
    len=atoi(argv[2]);
    pivot_strat = atoi(argv[3]);

    //only root 0 malloc memory for data, save the memory
    if (rank==0){
    	data = gen_num_vector(seq,len);

        // Prints the sequence.
        /*
    	printf("\n");
    	for (int i=0; i<len; i++) {
    		printf("%lf ", data[i]);
    	}
    	printf("\n");*/

    }

    /*** 1. Divide and distribute the data into p equal parts, one per process ***/
    chunk = len / size;
    remainder = len%size;
    if (rank==0){
    	if (remainder !=0){
    		printf("Error:Remainder is not 0.\nAccording to step1:Divide and distribute the data into p equal parts, one per process\nlength needs to be divisible by size\nOr wait us implement processing the remainder in the future\n");
    		return -1;
    	}
    }
    /* Prints the values that we send */
    /*
    if (rank == 0){

        for (int i =0; i < len; ++i)
            printf("Element %d: [%f]\n",i, data[i]);
    }
    */

    rcv_buffer = (double*)malloc(chunk*sizeof(double));

    MPI_Barrier(MPI_COMM_WORLD);
    parr_runtime = MPI_Wtime();

    MPI_Scatter(data, chunk, MPI_DOUBLE, rcv_buffer, chunk, MPI_DOUBLE,0,MPI_COMM_WORLD);

    /* Prints the values each PROC receives */
    //for (int i =0; i < chunk; ++i)
    //    printf("RANK[%d]\tdata[%d]=[%f]\n",rank, i, rcv_buffer[i]);

    /*** 2. Sort the data locally for each process ***/
    /* Use cmpfunc_desc for initial local sort if seq == -1 */
    if (seq==-1)
        qsort(rcv_buffer, chunk, sizeof(double), cmpfunc_desc); /* qsort() from stdlib does that for us */
    else
        qsort(rcv_buffer, chunk, sizeof(double), cmpfunc);

    /* Check that qsort() is successfull */
    /*Prints the values each PROC after local sort*/
    /*
    for (int i =0; i < chunk; ++i)
        printf("RANK[%d]\tElement %d: [%f]\n",rank, i, rcv_buffer[i]);
    */

    /*** 3. Perform global sort  ***/

    //Quicksort is recursive so it needs to be a function
    //parallel_quicksort(MPI_Comm comm, int rank, int size, double* data, int chunk);
    double* finaldata;
    if (size > 1){
	    loc_runtime = MPI_Wtime();
    	rcv_buffer = parallel_quicksort(MPI_COMM_WORLD,rank, size, rcv_buffer,chunk);
	    loc_runtime = MPI_Wtime() - loc_runtime;
    	//printf("\ncurrentsize:%d\n",currentsize);
        if( rank != 0 ){
            // All processes send to Zero.
            //MPI_Send( &currentsize, 1, MPI_INT, 0, rank, MPI_COMM_WORLD );
            MPI_Send( rcv_buffer, currentsize, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD );
        }else{
            // Zero process gathers data.
            int currentlocation = 0;
            finaldata=(double *)malloc(len*sizeof(double));
            memcpy(finaldata+currentlocation,rcv_buffer,currentsize*sizeof(double));
            currentlocation = currentlocation+currentsize;
            //printf("\n currentlocation:%d \n",currentlocation);
            double* buffer;
            //double* temp;
            for( int r = 1; r < size; ++r )
            {
                int bufferSize = 0;
                MPI_Probe(r, r, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, MPI_DOUBLE, &bufferSize);
                //MPI_Recv( &bufferSize, 1, MPI_INT, r, r, MPI_COMM_WORLD, NULL );
                //printf("\nbuffersize:%d\n",bufferSize);
                buffer = (double*)malloc(bufferSize*sizeof(double));
                MPI_Recv(buffer, bufferSize, MPI_DOUBLE, r, r, MPI_COMM_WORLD, NULL );
                memcpy(finaldata+currentlocation,buffer,bufferSize*sizeof(double));
                currentlocation = currentlocation+bufferSize;
                /*
    		    printf("\n");
    		    for (int i=0; i<bufferSize; i++) {
    		    	printf("buffer[%d] = %f, ", i, buffer[i]);
   		        }
    		    printf("\n");
    		*/   
                free(buffer);
            }
        }
    }



    MPI_Barrier(MPI_COMM_WORLD);
    parr_runtime = MPI_Wtime() - parr_runtime;
    double max_parr_runtime, max_loc_runtime;
    MPI_Reduce(&parr_runtime, &max_parr_runtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&loc_runtime, &max_loc_runtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        printf("%f\n", max_parr_runtime);
        //printf("%f\n", max_loc_runtime);
    }
    if (rank==0){
        /*
    	printf("\n");
    	for (int i=0; i<len; i++) {
    		printf("%f ", finaldata[i]);
    	}
    	printf("\n");
        */
        // Check results
    	int OK=1;
    	for (int i=0; i<len-1; i++) {
        	if(finaldata[i] > finaldata[i+1]) {
            		printf("Wrong:finaldata[%d] = %f, finaldata[%d] = %f\n", i, finaldata[i], i+1, finaldata[i+1]);
            		OK=0;
        	}
    	}
    	if (OK) printf("Data sorted correctly!\n");
    }
    //check_pivot(data,len,l_pivot,rank);
    if (rank==0) {
        free(data);
    }
    free(rcv_buffer);
    MPI_Finalize();
    return 0;
}

double* parallel_quicksort(MPI_Comm group, int rank, int size, double* data, int chunk){
    MPI_Comm_size( group, &size);
    MPI_Comm_rank( group, &rank );
    //Repeat recursively until each group consists of one single process.
    if (size < 2){
    	return (data);
    }
    /*
3.4 Merge the two sets of numbers in each process into one sorted list*/
    //3.1 Select pivot element within each process set
    double localpivot = pivot(pivot_strat, chunk, data, rank, size, group);

    //3.2 Locally in each process, divide the data into two sets according to the pivot (smaller or larger)
    double *low=(double *)malloc(chunk*sizeof(double));
    double *high=(double *)malloc(chunk*sizeof(double));
    int l ,h;
    l = 0;
    h = 0;
    for (int i = 0; i < chunk; i++){
    	if (data[i] <= localpivot){
    		low[l] = data[i];
    		l = l+1;
    	}else{
    		high[h] = data[i];
    		h = h+1;
    	}
    }
    //printf("\nl:%d,h:%d\n",l,h);
    //realloc could be used here

    //3.3 Split the processes into two groups and exchange data pairwise betweenthem so that all processes in one group get data less than the pivot and the others get data larger than the pivot.
    //use 8 processes as example
    //group1:0,1,2,3 gourp2:4,5,6,7
    int subsize = size/2;
    int pairrankh,pairrankl;
    int otherl,otherh;
    //int otherlt,otherht;
    if (rank < subsize){
    	pairrankh = rank + subsize;
    	otherl = 0;
    }else{
    	pairrankl = rank - subsize;
    	otherh = 0;
    }

    int tag = size;
    /*
    //Communicate size
    //This could be written using 1 MPI_Sendrecv, I split them to make the logic clear.
    if (rank < subsize){
    	//low process sends its highsize to pairhigh and receive pairhigh's lowsize from pairhigh
    	MPI_Sendrecv(&h, 1, MPI_INT, pairrankh, tag, &otherlt, 1, MPI_INT, pairrankh, tag+1, group, &status);

    }else{
    	//high process sends its lowsize to pairlow and receive pairlow's highsize from pairlow
    	MPI_Sendrecv(&l, 1, MPI_INT, pairrankl, tag+1, &otherht, 1, MPI_INT, pairrankl, tag, group, &status);


    }*/
	
    //Communicate high/low data
    MPI_Request request1;
    MPI_Request request2;
    MPI_Status status1;
    MPI_Status status2;
    double *otherlow;
    double *otherhigh;
    if (rank < subsize){
    	//low process sends its highdata to pairhigh and receive pairhigh's lowdata from pairhigh
    	
    	MPI_Isend(high, h, MPI_DOUBLE, pairrankh, tag+1, group,&request1);
    	MPI_Probe(pairrankh, tag, group, &status);
    	MPI_Get_count(&status, MPI_DOUBLE, &otherl);
    	//printf("\notherl:%d,otherlt:%d\n",otherl,otherlt);
    	otherlow=(double *)malloc(otherl*sizeof(double));
    	MPI_Irecv(otherlow, otherl, MPI_DOUBLE, pairrankh, tag, group, &request2);
    	MPI_Wait(&request1, &status1);
    	MPI_Wait(&request2, &status2);
    }else{
    	//high process sends its lowdata to pairlow and receive pairlow's highdata from pairlow
    	MPI_Isend(low, l, MPI_DOUBLE, pairrankl, tag, group,&request1);
    	MPI_Probe(pairrankl, tag+1, group, &status);
    	MPI_Get_count(&status, MPI_DOUBLE, &otherh);
    	//printf("\notherl:%d,otherlt:%d\n",otherh,otherht);
    	otherhigh=(double *)malloc(otherh*sizeof(double));
    	MPI_Irecv(otherhigh, otherh, MPI_DOUBLE, pairrankl, tag+1, group, &request2);
    	MPI_Wait(&request1, &status1);
    	MPI_Wait(&request2, &status2);	
    }

    //3.4 Merge the two sets of numbers in each process into one sorted list
    double* newdata;
    int newdatasize;
    free(data);
    MPI_Comm newgroup;
    int newsize,newrank;

    if (rank < subsize){
    	newdata = merge(low,l,otherlow,otherl);
    	newdatasize = l + otherl;
    	currentsize = newdatasize;
    	qsort(newdata, newdatasize, sizeof(double), cmpfunc);
    	/*
    	printf("\n");
    	for (int i=0; i<newdatasize; i++) {
    		printf("newdatah[%d] = %f, ", i, newdata[i]);
    	}
    	printf("\n");
    	*/
    	free(low);
    	free(otherlow);
    	// Split group into two. This is color 0.
        MPI_Comm_split( group, 0, 1, &newgroup );
        MPI_Comm_size( newgroup, &newsize );
        MPI_Comm_rank( newgroup, &newrank );
 	//Recursive call
 	return (parallel_quicksort(newgroup, newrank, newsize, newdata, newdatasize));
    }else{
    	newdata = merge(high,h,otherhigh,otherh);
    	newdatasize = h + otherh;
    	currentsize = newdatasize;
    	qsort(newdata, newdatasize, sizeof(double), cmpfunc);
    	/*
    	printf("\n");
    	for (int i=0; i<newdatasize; i++) {
    		printf("newdatah[%d] = %f, ", i, newdata[i]);
    	}
    	printf("\n");
    	*/
    	free(high);
    	free(otherhigh);
    	// Split group into two. This is color 1.
        MPI_Comm_split( group, 1, 1, &newgroup );
        MPI_Comm_size( newgroup, &newsize );
        MPI_Comm_rank( newgroup, &newrank );
        //Recursive call
        return (parallel_quicksort(newgroup, newrank, newsize, newdata, newdatasize));
    }

}//quicksort parallel

double pivot(int pivot_strat,int len, double *rcv_buffer, int rank, int size, MPI_Comm group){
    double l_pivot = 0;
    MPI_Comm_size( group, &size);
    MPI_Comm_rank( group, &rank );
    switch (pivot_strat){
        case STRAT_ONE : /* 1. Select the median in one processor in each group of processors*/
        {
            /* ROOT find the Global PIVOT and broadcast it to all other PROCS */
            if (rank==0){
                l_pivot = find_pivot(rcv_buffer, len);
                //l_pivot = find_mean(data, len);

                //for (int i =0; i < chunk; ++i)
                //    printf("RANK[%d]\tElement %d: [%f]\n",rank, i, data[i]);

                //printf("Global Pivot [%f]\n", l_pivot);

            }
            //MPI_Bcast ( &l_pivot, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast ( &l_pivot, 1, MPI_DOUBLE, 0, group);
            //printf("Global Pivot [%f]", l_pivot);
	    return l_pivot;
        }
        break;
        case STRAT_TWO : /* Select the median of all medians in each processor group */
        {
            double local_pivot = find_pivot(rcv_buffer, len); /* Each PROC find median of local buffer */
            double all_median[size]; /* # of median values as our #PROCS */
            //MPI_Allgather( &local_pivot, 1, MPI_DOUBLE, all_median, 1, MPI_DOUBLE, MPI_COMM_WORLD);
            MPI_Allgather( &local_pivot, 1, MPI_DOUBLE, all_median, 1, MPI_DOUBLE, group);
            l_pivot = find_pivot(all_median, size); /* median of all medians */
            //printf("Global Pivot [%f]\n", l_pivot);
	    return l_pivot;
        }
        break;
        case STRAT_THREE: /* Select the mean value of all medians in each processor group */
        {
            double local_pivot = find_pivot(rcv_buffer, len); /* Each PROC find median of local buffer */
            double all_median[size]; /* # of median values as our #PROCS */
            //MPI_Allgather( &local_pivot, 1, MPI_DOUBLE, all_median, 1, MPI_DOUBLE, MPI_COMM_WORLD);
            MPI_Allgather( &local_pivot, 1, MPI_DOUBLE, all_median, 1, MPI_DOUBLE, group);
            l_pivot = find_mean( all_median, size);
            //printf("Global Pivot [%f]\n", l_pivot);
	    return l_pivot;
        }
        break;
    }
    return 0;
}

double find_pivot(double* data, int len){
    /* If the number of elements is zero the pivot is data[0] */
    if (len==0)
        return 0;
    double pivot = 0;
    /* If number of elements is odd, pivot is the middle element */
    if (len % 2 != 0)
        pivot = data[len/2];
    /* Else it is the average of the two middle elements */
    else
        pivot =  ((double)data[(len / 2) -1] + (double)data[(len / 2)]) / 2.0;
    return pivot;
}

double find_mean(double *data, int len){
    double sum = 0;
    for (int i=0; i<len; ++i)
        sum += data[i];
    sum = sum/(double)len;
    return sum;
}

void check_pivot(double* data, int len, double pivot, int rank){
    int lower,upper;
    for (int i=0; i<len; i++){
        if (data[i] > pivot)
            upper++;
        else
            lower++;

    }
    if (rank == 0)
        printf("#lower=%d\t#upper=%d\n", lower, upper);
}

double* gen_num_vector(int seq, int len){

    double *data=(double *)malloc(len*sizeof(double));
    int i;
    if (seq==0) {

        // Uniform random numbers
        for (i=0;i<len;i++)
        data[i]=drand48();

    }

    else if (seq==1) {

        // Exponential distribution
        double lambda=10;
        for (i=0;i<len;i++)
        data[i]=-lambda*log(1-drand48());
    }

    else if (seq==2) {

        // Normal distribution
        double x,y;
        for (i=0;i<len;i++){
            x=drand48(); y=drand48();
            data[i]=sqrt(-2*log(x))*cos(2*PI*y);
        }
    }

    else if (seq==-1) {
	// Descending order
	for (i=0;i<len;i++) {
	    data[i] = len-i;
	}
    }

    return data;
}

int cmpfunc_desc(const void * a, const void * b){

    if (*(double*)a > *(double*)b)
        return -1;
    else if (*(double*)a < *(double*)b)
        return 1;
    else
        return 0;
}

int cmpfunc (const void * a, const void * b)
{
    if (*(double*)a > *(double*)b)
        return 1;
    else if (*(double*)a < *(double*)b)
        return -1;
    else
        return 0;
}

double *merge(double *v1, int n1, double *v2, int n2)
{
    int i,j,k;
    double *result;

    result = (double *)malloc((n1+n2)*sizeof(double));

    i=0; j=0; k=0;
    while(i<n1 && j<n2)
        if(v1[i]<v2[j])
        {
            result[k] = v1[i];
            i++; k++;
        }
        else
        {
            result[k] = v2[j];
            j++; k++;
        }
    if(i==n1)
        while(j<n2)
        {
            result[k] = v2[j];
            j++; k++;
        }
    else
        while(i<n1)
        {
            result[k] = v1[i];
            i++; k++;
        }
    return result;
}
