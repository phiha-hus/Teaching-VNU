CC = gcc
RM = rm -f

CFLAGS = -O2 -march=native -funroll-loops
LIBS = -lm
LDFLAGS = 

SOURCES_LIB :=	\
	basic.c \
	linalg.c \
	ode.c

HEADER_LIB := \
	$(SOURCES_LIB:.c=.h)

OBJECTS_LIB := \
	$(SOURCES_LIB:.c=.o)	

SOURCES_PROGRAMS :=	\
	exercise_mass_spring.c \
	exercise_gravitation.c

PROGRAMS := \
	$(SOURCES_PROGRAMS:.c=)

all: $(PROGRAMS)

$(OBJECTS_LIB): %.o: %.c
	@echo Compiling \"$<\"
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGRAMS): %: %.c $(OBJECTS_LIB)
	@echo Compiling and linking \"$<\"
	$(CC) $< $(HEADER_LIB) -o $@ $(LDFLAGS) $(LIBS) $(OBJECTS_LIB)

clean:
	$(RM) $(PROGRAMS) $(OBJECTS_LIB)