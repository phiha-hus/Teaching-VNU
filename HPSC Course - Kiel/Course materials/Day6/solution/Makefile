CC = gcc
RM = rm -f

CFLAGS = -O3 -march=native -funroll-loops -fopenmp
LIBS = -lm
LDFLAGS = -fopenmp

SOURCES_LIB :=	\
	basic.c \
	linalg.c \
	miniblas.c \
	gridfunc2d.c

HEADER_LIB := \
	$(SOURCES_LIB:.c=.h)

OBJECTS_LIB := \
	$(SOURCES_LIB:.c=.o)	

SOURCES_PROGRAMS :=	\
	exercise_omp_gridfunc2d.c \
	exercise_montecarlo.c \
	exercise_par_block_rl.c

PROGRAMS := \
	$(SOURCES_PROGRAMS:.c=)

all: $(PROGRAMS)

$(OBJECTS_LIB): %.o: %.c
	@echo Compiling \"$<\"
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGRAMS): %: %.c $(OBJECTS_LIB)
	@echo Compiling and linking \"$<\"
	$(CC) $(CFLAGS) $< $(HEADER_LIB) -o $@ $(LDFLAGS) $(OBJECTS_LIB) $(LIBS)

clean:
	$(RM) $(PROGRAMS) $(OBJECTS_LIB)
