# Readme
## Advice:
According to Jarmo's feedback:  
### Main: 
1.When you collect the data on sender side send only the data and not the size. On rank 0 use MPI_Probe, MPI_Get_size and MPI_Recv.***（Finished by Junjie 2021/May/31）***     
2.Then DON'T use merge, just copy the data into the final array.***（Finished by Junjie 2021/May/31）***   
### parallel_quicksort: 
1.You don't have to allocate low, high array. just find the split point by searching either from the middle with linear search or by using binary search.***（Finished by Junjie 2021/June/01）***     
2.don't exchange size information, use Probe and get_size functions as above.***（Finished by Junjie 2021/May/31）***   
3.Use MPI_Isend to send the data and MPI_Recv to receive. After exchange you should use the merge function.***（Finished by Junjie 2021/May/31）***   
4.Why do you qsort on lines 301 and 323, is this needed, how does your data look after the merge? It seems that you miss the point with this algorithm when you start sorting in each step.***（comment out the 2 qsort! really no need to do that because merge has made the data sorted! Finished by Junjie 2021/June/13）***       
5.When you collect the data why do you need an extra buffer, why not receiving directly into the final array?***（Modify the code and store the sorted data from different processors directly in the finaldata! Finished by Junjie 2021/June/13）***    

# *Update:2021/June/13 Junjie Chu*  
To complie:  
```
mpicc -std=c99 -g -O3 -o quicksort_revised_0 quicksort_revised_0.c -lm
```
## Results:
### 2021/6/01 Junjie
![image](https://user-images.githubusercontent.com/65893273/121814468-d928d700-cca3-11eb-981e-068a1b59027a.png)    
speed up the code!!!  



# *Update:2021/June/01 Junjie Chu*  
To complie:  
```
mpicc -std=c99 -g -O3 -o quicksort_revised quicksort_revised.c -lm
```
## Results:
### 2021/6/01 Junjie
![image](https://user-images.githubusercontent.com/65893273/120345782-5e6bcd80-c32d-11eb-8d60-455427a874e1.png)  

# *Update:2021/May/31 Junjie Chu*  
To complie:  
```
mpicc -std=c99 -g -O3 -o quickrevised quickrevised.c -lm
```

## Results:
### 2021/5/31 Junjie
![image](https://user-images.githubusercontent.com/65893273/120116501-0e6cf980-c1bb-11eb-9a48-b47b10d7f0bb.png)  
