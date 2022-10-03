/* File: mpi_vector_ops.c
 * COMP 137-1 Spring 2020
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h> 
char* infile = NULL;
char* outfile = NULL;

int readInputFile(char* filename, long* n_p, double* x_p, double** A_p, double** B_p)
{
    long i;
    FILE* fp = fopen(filename, "r");
    if (fp == NULL) return 0;
    fscanf(fp, "%ld\n", n_p);
    fscanf(fp, "%lf\n", x_p);
    *A_p = malloc(*n_p*sizeof(double));
    *B_p = malloc(*n_p*sizeof(double));
    for (i=0; i<*n_p; i++) fscanf(fp, "%lf\n", (*A_p)+i);
    for (i=0; i<*n_p; i++) fscanf(fp, "%lf\n", (*B_p)+i);
    return 1;
}

int writeOutputFile(char* filename, long n, double y, double* C, double* D)
{
    long i;
    FILE* fp = fopen(filename, "w");
    if (fp == NULL) return 0;
    fprintf(fp, "%ld\n", n);
    fprintf(fp, "%lf\n", y);
    for (i=0; i<n; i++) fprintf(fp, "%lf\n", C[i]);
    for (i=0; i<n; i++) fprintf(fp, "%lf\n", D[i]);
    return 1;
}

void serialSolution(long n, double x, double* A, double* B,
                    double* y, double* C, double* D)
{
    long i;
    double dp = 0.0;
    for (i=0; i<n; i++)
    {
        C[i] = x*A[i];
        D[i] = x*B[i];
        dp += A[i]*B[i];
    }
    *y = dp;
}

int main(int argc, char* argv[])
{
    long    n=0;     /* size of input arrays */
    double  x;       /* input scalar */
    double* A;       /* input vector */
    double* B;       /* input vector */
    double* C;       /* output vector xA */
    double* D;       /* output vector xB */
    double  y;       /* output scalar A dot B */
    int my_rank; 
    int comm_sz; 
    double* totalc; 
    double* totald; 
    double total = 0; 

    MPI_Init(NULL, NULL); 
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz); 
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 
 
    if (my_rank == 0) {
              
    /* read input data */
    if (argc<3)
    {
        n = -1;
        fprintf(stderr, "Command line arguments are required.\n");
        fprintf(stderr, "argv[1] = name of input file\n");
        fprintf(stderr, "argv[2] = name of input file\n");
    }
    else
    {
        infile = argv[1];
        outfile = argv[2];
        if (!readInputFile(infile, &n, &x, &A, &B))
        {
            fprintf(stderr, "Error opening input files. Aborting.\n");
            n = -1;
        }
    }
}
    if (n < 0)
    {
        fprintf(stderr, "Aborting task due to input errors.\n");
        exit(1);
    }    
    
    MPI_Bcast(&x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    long section = n / comm_sz; 
    
    if(n % comm_sz !=0){
    printf("processes and size are not correct");
    exit(1);
    }

    double* APartition = malloc(section * sizeof(double));    
    double* BPartition = malloc(section * sizeof(double));
    C = malloc((section) * sizeof(double));   
    D = malloc((section) * sizeof(double)); 
    
    
    double sub_total = 0; 
     
     MPI_Scatter(A, section, MPI_DOUBLE, APartition, section, MPI_DOUBLE, 0, MPI_COMM_WORLD);     
    MPI_Scatter(B, section, MPI_DOUBLE, BPartition, section, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    serialSolution(section, x, APartition, BPartition, &sub_total, C, D);  
 MPI_Reduce(&sub_total, &total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    printf("total %d", total); 
    /* for (i = 0; i < n; i++) {
       C[i] = x*CP[i];
       D[i] = x*DP[i]; 
       dp += A[i]*B[i]; 
    } */ 
    totalc = malloc(n * sizeof(double)); 
    totald = malloc(n * sizeof(double));
    
    if (my_rank == 0) {
	MPI_Gather(C, section, MPI_DOUBLE, totalc, section, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather(D, section, MPI_DOUBLE, totald, section, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    
    else {
    	MPI_Gather(C, section, MPI_DOUBLE, totalc, section, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    	MPI_Gather(D, section, MPI_DOUBLE, totald, section, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } 

    if (my_rank == 0) {
      

    if (!writeOutputFile(outfile, n, total, totalc, totald))
    {
        fprintf(stderr, "Error opening output file. Aborting.\n");
        exit(1);
    }
}
    /* free all dynamic memory allocation */
    //free(A);
    //free(B);
    free(C);
    free(D);
    MPI_Finalize(); 
    return 0;
}

