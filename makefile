mpi_vector_ops: mpi_vector_ops.c
	mpicc -g mpi_vector_ops.c -o mpi_vector_ops

clean:
	rm mpi_vector_ops
