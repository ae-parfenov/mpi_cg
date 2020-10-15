CPP = mpicxx
FLAGS = -std=gnu++11

all: mpi_cg

mpi_cg: main.o global_entities.o spaces_connection.o local_entities.o linalg_operations.o
	$(CPP) main.o global_entities.o spaces_connection.o local_entities.o linalg_operations.o -o mpi_cg

main.o: main.cpp
	$(CPP) -c main.cpp -o main.o $(CXXFLAGS)

global_entities.o: global_entities.cpp
	$(CPP) -c global_entities.cpp -o global_entities.o $(FLAGS)

spaces_connection.o: spaces_connection.cpp
	$(CPP) -c spaces_connection.cpp -o spaces_connection.o $(FLAGS)

local_entities.o: local_entities.cpp
	$(CPP) -c local_entities.cpp -o local_entities.o $(FLAGS)

linalg_operations.o: linalg_operations.cpp
	$(CPP) -c linalg_operations.cpp -o linalg_operations.o $(FLAGS)

shared_operations.o: shared_operations.cpp
	$(CPP) -c shared_operations.cpp -o shared_operations.o $(FLAGS)

clean:
	rm -rf mpi_cg main.o global_entities.o spaces_connection.o local_entities.o linalg_operations.o
