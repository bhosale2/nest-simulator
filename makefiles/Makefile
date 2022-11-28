hn = $(shell hostname)
username = $(shell whoami)

# Set this one to your own custom path
CC=/usr/local/Cellar/gcc/12.1.0/bin/g++-12

CPPSETTINGS += -std=c++11 
CPPFLAGS+= -O3 -Wno-all -Wno-deprecated -w -Wno-reorder -DNDEBUG -march=native

NEST_OBJS = \
	main.o \
	GeometryFunctions.o \
	Polymer.o \
	PolymerIntegrator.o \
	PositionVerlet2nd.o \
	Rod.o \
	RodInitialConfigurations.o \
	Vector3.o \
	Matrix3.o \
	SpeedFunctions.o \
	MathFunctions.o \
	Nest.o \
	Mindlin_friction.o \
	Collision_test.o \

VPATH := ../source/
.DEFAULT_GOAL := nest

#Linking stage
nest: $(NEST_OBJS)
	$(CC) $(CPPSETTINGS) $^ $(LIB) -o  $@
	@echo done

#Compiling stage	
%.o:%.cpp
	$(CC) $(CPPSETTINGS) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o
	rm -f nest
