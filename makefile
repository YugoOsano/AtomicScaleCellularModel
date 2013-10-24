
CC = g++
CFLAGS = -c -Wall -Wno-deprecated -ansi -O2

OFILES = monte_carlo.o random_number.o  rotating_matrix.o \
	shape_trim.o   shape_remove_isolation.o  shape_trim2.o\
	shape_stickingCl.o \
	particle.o  neutral_particle.o  ion_particle.o \
	integration_angle.o  \
	fileio_particle.o  iedf_to_distribution.o \
	common_utils.o

all: monte_carlo

monte_carlo: $(OFILES)
	$(CC) $(OFILES) -o monte_carlo

monte_carlo.o: monte_carlo.cc atom_struct.h header_main.h\
	random_number.h \
	integration_angle.cc rotating_matrix.h shape_trim.h \
	particle.h  neutral_particle.h   ion_particle.h\
	integration_angle.h  fileio_particle.h  common_utils.h
	$(CC) $(CFLAGS) monte_carlo.cc 

rotating_matrix.o: rotating_matrix.cc rotating_matrix.h \
	atom_struct.h 
	$(CC) $(CFLAGS) rotating_matrix.cc

shape_trim.o: shape_trim.cc shape_trim.h  atom_struct.h
	$(CC) $(CFLAGS) shape_trim.cc

shape_trim2.o: shape_trim2.cc shape_trim.h  atom_struct.h \
	fileio_particle.h
	$(CC) $(CFLAGS) shape_trim2.cc

shape_stickingCl.o: shape_stickingCl.cc  shape_trim.h
	$(CC) $(CFLAGS) shape_stickingCl.cc 

shape_remove_isolation.o: shape_remove_isolation.cc \
	shape_trim.h  atom_struct.h
	$(CC) $(CFLAGS)  shape_remove_isolation.cc

integration_angle.o: integration_angle.cc integration_angle.h
	$(CC) $(CFLAGS)  integration_angle.cc

particle.o: particle.cc  particle.h  rotating_matrix.h \
	integration_angle.h  random_number.h  atom_struct.h \
	iedf_to_distribution.h  fileio_particle.h  common_utils.h
	$(CC) $(CFLAGS) particle.cc  

neutral_particle.o: neutral_particle.h neutral_particle.cc\
	particle.h   atom_struct.h  \
	shape_trim.h
	$(CC) $(CFLAGS) neutral_particle.cc  


ion_particle.o:   ion_particle.cc   ion_particle.h  \
	particle.h  shape_trim.h \
	atom_struct.h
	$(CC) $(CFLAGS)    ion_particle.cc

random_number.o: random_number.cc random_number.h
	$(CC) $(CFLAGS) random_number.cc


fileio_particle.o: fileio_particle.h fileio_particle.cc
	$(CC) $(CFLAGS)  fileio_particle.cc

iedf_to_distribution.o: iedf_to_distribution.h iedf_to_distribution.cc
	$(CC) $(CFLAGS)  iedf_to_distribution.cc

common_utils.o: common_utils.cc  common_utils.h\
	atom_struct.h  fileio_particle.h   
	$(CC) $(CFLAGS)  common_utils.cc

clean:
	rm -f $(OFILES)
