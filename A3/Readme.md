# Readme
*Update:2021/May/31 Junjie Chu*

According to Jarmo's feedback:  
## Main: 
1.When you collect the data on sender side send only the data and not the size. On rank 0 use MPI_Probe, MPI_Get_size and MPI_Recv.***（Finished by Junjie）***     
2.Then DON'T use merge, just copy the data into the final array.***（Finished by Junjie）***   
## parallel_quicksort: 
1.You don't have to allocate low, high array. just find the split point by searching either from the middle with linear search or by using binary search.   
2.don't exchange size information, use Probe and get_size functions as above.***（Finished by Junjie）***   
3.Use MPI_Isend to send the data and MPI_Recv to receive. After exchange you should use the merge function.***（Finished by Junjie）***   
