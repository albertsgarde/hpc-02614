# Makefile
#
TARGET_J  = poisson_j		# Jacobi
TARGET_GS = poisson_gs		# Gauss-Seidel

SOURCES	= main.c print.c alloc3d.c
OBJECTS	= print.o alloc3d.o frobenius_norm.o init_u.o init_f.o utils.o test_solution.o
MAIN_J	= main_j.o
MAIN_GS = main_gs.o
OBJS_J	= jacobi_ser.o jacobi_par1.o jacobi_par2.o $(MAIN_J)
OBJS_GS	= gauss_seidel_ser.o gauss_seidel_par1.o gauss_seidel_par2.o $(MAIN_GS)

# options and settings for the GCC compilers
#
CC	= gcc
DEFS	= 
OPT	= -g -Ofast
IPO	= 
ISA	= 
CHIP	= 
ARCH	= 
PARA	= 
CFLAGS	= $(DEFS) $(ARCH) $(OPT) $(ISA) $(CHIP) $(IPO) $(PARA) $(XOPTS) "-std=gnu11" -fopenmp 
LDFLAGS = -lm 

all: $(TARGET_J) $(TARGET_GS) 

$(TARGET_J): $(OBJECTS) $(OBJS_J)
	$(CC) -o $@ $(CFLAGS) $(OBJS_J) $(OBJECTS) $(LDFLAGS)

$(TARGET_GS): $(OBJECTS) $(OBJS_GS)
	$(CC) -o $@ $(CFLAGS) $(OBJS_GS) $(OBJECTS) $(LDFLAGS)

$(MAIN_J):
	$(CC) -o $@ -D_JACOBI $(CFLAGS) -c main.c 

$(MAIN_GS):
	$(CC) -o $@ -D_GAUSS_SEIDEL $(CFLAGS) -c main.c 

clean:
	@/bin/rm -f core *.o *~

realclean: clean
	@/bin/rm -f $(TARGET_J)  $(TARGET_GS)

# DO NOT DELETE

main_j.o: main.c utils.h test_solution.h init_u.h init_f.h print.h jacobi_ser.h jacobi_par1.h jacobi_par2.h
main_gs.o: main.c utils.h test_solution.h init_u.h init_f.h print.h gauss_seidel_ser.h gauss_seidel_par1.h gauss_seidel_par2.h
print.o: print.h
