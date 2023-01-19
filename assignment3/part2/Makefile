TARGET= poisson
OBJS	= warm_up.o \
        utils.o \
        frobenius_norm.o \
        alloc3d.o \
        d_alloc3d.o \
        init_f.o \
        init_u.o \
        test_solution.o \
        print.o \
        jacobi_seq.o \
        jacobi_par.o \
        jacobi_gpu_map.o \
        jacobi_gpu_mcp.o \
        jacobi_gpu_async.o \
        main.o \

CC	= nvc
CXX	= nvc++

OPT	= -g -fast -Msafeptr -Minfo -acc -mp=gpu -gpu=pinned -gpu=cc80 -gpu=lineinfo -cuda -mp=noautopar
ISA	= 
PARA	= --diag_suppress deprecated_string_conv
INC   = -I$(CUDA_ROOT)/include -I/appl/nvhpc/2022_2211/Linux_x86_64/22.11/examples/OpenMP/SDK/include
LIBS	= -lcuda

CFLAGS = $(OPT) $(INC) $(ISA) $(PARA) $(XOPT)
CXXFLAGS = $(OPT) $(INC) $(ISA) $(PARA) $(XOPT)

all: $(TARGET)

$(TARGET): $(OBJS) 
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIBS)
clean:
	@/bin/rm -f *.o $(TARGET)
