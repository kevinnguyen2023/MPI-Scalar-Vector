# MPI Scalar Vector 
Developed an MPI program that computes two scalar vector products and one dot product in C. Implemented functions to read the input/output files along with makefile that properly distributes the data. Used MPI functions that distributes input data from the size of the arrays and computes the product of scalar vectors. 

# Program Requirements
- The program distributes the work evenly among p processes. 
- The program assumes the size of the arrays is divisible by the number of processes, n mod p = 0
- Process 0 reads the input and output file meaning it must distribute the input data to other p-1 processes. 
