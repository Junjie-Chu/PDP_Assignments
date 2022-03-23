There are 5 files(except the report) in the folder: topic2.c,prop.h,prop.c,Makefile,simulatemalaria and a folder including the results when N= 1000000,1100000,1200000
# Example:   
## Step 1: make all  
You will see a program named as simulatemalaria.  

## Step 2: mpirun -np 4 ./simulatemalaria 10000 result4*10000    
There are 2 arguments in the command: '10000' is the local data size of each processor and 'result4*10000' is the name of output file.  
In the case, we use 4 processors and each processor does 10000 experiments. So in total, 4*10000 experiments are simulated.  
If the total size is no larger than 400, a simple histogram is printed.  

## Step 3:
You will see an output file named as result4*10000.

The output file contains two lines: the 1st line is 21 integers indicating the range of histogram with even intervals seperated by ","; the 2nd line is 20 integers indicating the number of susceptible cases that belongs to each intervals respectively. 

*Author: Junjie Chu*
