
CXX = g++
CXXFLAGS = -Wall -g -I.

HDF5_ROOT_DIR = /usr/local/Cellar/hdf5/1.10.5
CCFLAGS_HDF5 = -I$(HDF5_ROOT_DIR)/include/
CLDFLAGS_HDF5 = -L$(HDF5_ROOT_DIR)/lib -lhdf5_hl -lhdf5

DEPS = grid.hpp data_column.hpp algorithms.hpp \
       right_hand_side.hpp time_stepping.hpp

DEPS += enumerators.hpp make_arrays.hpp \
	    reconstruction.hpp riemann_solvers.hpp \
	    tools.hpp parameters.hpp

OBJ = grid.o data_column.o algorithms.o \
      right_hand_side.o time_stepping.o main.o

mars: $(OBJ)
	$(CXX) $(CXXFLAGS) $(CCFLAGS_HDF5) $^ -o $@ $(CLDFLAGS_HDF5)

%.o: %.cpp $(DEPS)
	$(CXX) $(CCFLAGS_HDF5) -c -o $@ $<

clean:
	rm -f *.o
