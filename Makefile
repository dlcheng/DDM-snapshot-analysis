#------------------------------------------  Always recommended
OPT += -DNO_AHF_PROFILE
OPT += -DID_START=0                          # the particle ID starts from 0
#------------------------------------------- Decay dark matter staff
#OPT += -DDECAY_DARK_MATTER
#------------------------------------------- Used for debug
#OPT += -DDEBUG
#OPT += -DVERBOSE
#------------------------------------------  Select target computer
#SYSTYPE="dlcheng"
#SYSTYPE="ITSC"
SYSTYPE="Mac"
#------------------------------------------  Adjust settings for target computer

ifeq ($(SYSTYPE),"dlcheng")
CC       =  gcc   
OPTIMIZE =  -O3 -Wall
GSL_LIBS=   -L/Users/dalongcheng/Install/gsl-1.16/lib
GSL_INCL =  -I/Users/dalongcheng/Install/gsl-1.16/include
endif

ifeq ($(SYSTYPE),"ITSC")
CC       =  gcc   
OPTIMIZE =  -O3 -Wall
GSL_LIBS =  -L/users/s0902248/Lib/gsl-1.9/lib  -Wl,"-R /users/s0902248/Lib/gsl-1.9/lib"
GSL_INCL =  -I/users/s0902248/Lib/gsl-1.9/include
endif

ifeq ($(SYSTYPE),"Mac")
CC       =  gcc   
OPTIMIZE =  -O3 -Wall
GSL_LIBS=   -L/usr/local/Cellar/gsl/1.16/lib 
GSL_INCL =  -I/usr/local/Cellar/gsl/1.16/include
endif

OPTIONS  =  $(OPTIMIZE) $(OPT) 

EXEC     =  Halo_info

OBJS     =  allvars.o construct_halos.o file_close.o gsl_solver.o \
            halo_profile.o init_all.o load_ahf_file.o \
            load_gadget.o main.o \
            memory_control.o  warn_error.o \
            write_file.o

INCL     =  allvars.h proto.h define.h Makefile

CFLAGS   =  $(OPTIONS) $(GSL_INCL)

LIBS     =  $(GSL_LIBS) -lgsl -lgslcblas -lm 

$(EXEC): $(OBJS) 
	$(CC) $(OBJS) $(LIBS)  -o  $(EXEC)  

$(OBJS): $(INCL) 

clean:
	rm -f $(OBJS) $(EXEC)