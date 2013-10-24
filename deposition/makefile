
CC = g++
#CFLAGS = -c -Wno-deprecated -ansi -O2
CFLAGS = -c -Wall -Wno-deprecated -ansi -O2 -I../

OFILES = main_deposition.o random_number.o  rotating_matrix.o \
	shape_trim.o    shape_trim2.o \
	particle.o  neutral_particle.o etch_products.o ion_particle.o\
	integration_angle.o  \
	fileio_particle.o  iedf_to_distribution.o \
	common_utils.o 

all : main_deposition

main_deposition: $(OFILES)
	$(CC) $(OFILES) -o      main_deposition

#monte_carlo.o: ../monte_carlo.cc ../atom_struct.h ../header_main.h\
#	header_deposition.h   ../random_number.h \
#	../integration_angle.cc ../rotating_matrix.h ../shape_trim.h \
#	../particle.h  ../neutral_particle.h \
#	../integration_angle.h  ../fileio_particle.h
#	$(CC) $(CFLAGS) ../monte_carlo.cc 

main_deposition.o: main_deposition.cc  ../monte_carlo.cc \
	included_to_main1.cc   included_to_main2.cc \
	../atom_struct.h    header_deposition.h  \
	../header_main.h   ../fileio_particle.h  \
	../particle.h   ../neutral_particle.h 	etch_products.h\
	../ion_particle.h \
	../shape_trim.h   ../common_utils.h
	$(CC) $(CFLAGS) main_deposition.cc

rotating_matrix.o: ../rotating_matrix.cc ../rotating_matrix.h \
	../atom_struct.h 
	$(CC) $(CFLAGS) ../rotating_matrix.cc

shape_trim.o: ../shape_trim.cc ../shape_trim.h  ../atom_struct.h
	$(CC) $(CFLAGS) ../shape_trim.cc

shape_trim2.o: ../shape_trim2.cc ../shape_trim.h  ../atom_struct.h \
	../fileio_particle.h
	$(CC) $(CFLAGS) ../shape_trim2.cc

#shape_stickingCl.o: shape_stickingCl.cc  shape_trim.h
#	$(CC) $(CFLAGS) shape_stickingCl.cc 

#shape_remove_isolation.o: shape_remove_isolation.cc \
#	shape_trim.h  atom_struct.h
#	$(CC) $(CFLAGS)  shape_remove_isolation.cc

integration_angle.o: ../integration_angle.cc ../integration_angle.h
	$(CC) $(CFLAGS)  ../integration_angle.cc

particle.o: ../particle.cc  ../particle.h  ../rotating_matrix.h \
	../integration_angle.h  ../random_number.h  ../atom_struct.h \
	../iedf_to_distribution.h  ../fileio_particle.h  \
	../common_utils.h
	$(CC) $(CFLAGS) ../particle.cc  

neutral_particle.o: ../neutral_particle.h ../neutral_particle.cc\
	../particle.h   ../atom_struct.h  \
	../shape_trim.h
	$(CC) $(CFLAGS) ../neutral_particle.cc  

etch_products.o:  etch_products.cc  etch_products.h \
	../neutral_particle.h ../neutral_particle.cc\
	../particle.h   ../atom_struct.h ../header_main.h
	$(CC) $(CFLAGS)  etch_products.cc

ion_particle.o:  ../ion_particle.cc  ../ion_particle.h \
	../iedf_to_distribution.h  ../random_number.h  ../particle.h  \
	../shape_trim.h  ../atom_struct.h
	$(CC) $(CFLAGS)   ../ion_particle.cc

random_number.o: ../random_number.cc ../random_number.h
	$(CC) $(CFLAGS) ../random_number.cc


fileio_particle.o: ../fileio_particle.h ../fileio_particle.cc
	$(CC) $(CFLAGS)  ../fileio_particle.cc

iedf_to_distribution.o: ../iedf_to_distribution.h ../iedf_to_distribution.cc
	$(CC) $(CFLAGS)  ../iedf_to_distribution.cc

common_utils.o: ../common_utils.cc ../common_utils.h\
	../atom_struct.h   ../fileio_particle.h   
	$(CC) $(CFLAGS)  ../common_utils.cc

clean:
	rm -f $(OFILES)
