#
CC=gcc
COPTS= -O2 -Wall


libmesh6.o :	libmesh6.c libmesh6.h 
	$(CC) -c $(COPTS)  -I. libmesh6.c  

mesh.o :	mesh.c mesh.h libmesh6.h 
	$(CC) -c $(COPTS)  -I. mesh.c

eigen.o :	eigen.c mesh.h 
	$(CC) -c $(COPTS)  -I. eigen.c


# -------------------------- TP1 ---------------
mesh:	mesh.o main_mesh.o eigen.o libmesh6.o  
	$(CC) $(COPTS) -o mesh mesh.o  main_mesh.o eigen.o libmesh6.o -lm

main_mesh.o :	main_mesh.c mesh.h
	$(CC) -c $(COPTS)  -I. main_mesh.c  

# ------------------------- TP1 ex3 ----------------

domain: mesh.o main_domain.o eigen.o libmesh6.o  
	$(CC) $(COPTS) -o domain mesh.o  main_domain.o eigen.o libmesh6.o -lm

main_domain.o :	main_domain.c mesh.h
	$(CC) -c $(COPTS)  -I. main_domain.c  
# ------------------------- TP2 delauney ---------------- 
cavity.o :	main_cavity.c mesh.h
	$(CC) -c $(COPTS)  -I. main_cavity.c 

cavity:	mesh.o cavity.o eigen.o libmesh6.o  
	$(CC) $(COPTS) -o cavity mesh.o  main_cavity.o eigen.o libmesh6.o -lm 

# ------------------------- TP2 compression -------------

compression: mesh.o main_compression.o eigen.o libmesh6.o  
	$(CC) $(COPTS) -o compression mesh.o  main_compression.o eigen.o libmesh6.o -lm 
	
main_compression.o :	main_compression.c mesh.h
	$(CC) -c $(COPTS)  -I. main_compression.c 

# ------------------------------------------------------

clean :
	-rm mesh.o main_mesh.o eigen.o libmesh6.o


