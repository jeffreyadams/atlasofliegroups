# sources contains a list of the source files (i.e., the .cpp files)
sources := $(patsubst %.cpp,%.cpp,$(wildcard sources/*/*.cpp))
# there is one .o file for each .cpp file
objects := $(patsubst %.cpp,%.o,$(wildcard sources/*/*.cpp))
# this is used to automatically generate the dependencies
dependencies := $(patsubst %.cpp,%.d,$(wildcard sources/*/*.cpp))

includedirs := $(addprefix -I,$(wildcard sources/*))

pflags = -c $(includedirs) -pg -O -DNREADLINE
oflags = -c $(includedirs) -O -Wall -DNDEBUG
gflags = -c $(includedirs) -g

cflags := $(oflags) # the default setting

ifeq ($(debug),true)
	cflags := $(gflags)
endif

ifeq ($(optimize),true)
	cflags := $(oflags)
endif

ifeq ($(profile),true)
	cflags := $(pflags)
endif

ifeq ($(readline),false)
	cflags := $(cflags) -DNREADLINE
else
	rlincludes := -lreadline -lcurses
#use this for readline on the Mac. 
#       rlincludes := -lreadline.5.1 -lcurses
endif

ifeq ($(verbose),true)
	cflags := $(cflags) -DVERBOSE
endif

cc = g++ # the default compiler

# give compiler=icc argument to make to use the intel compiler
ifdef compiler
	cc = $(compiler)
endif

all: atlas # clean

atlas: $(objects)
ifeq ($(profile),true)
	$(cc) -pg -o atlas $(objects)
else
	$(cc) -o atlas $(objects) $(rlincludes)
endif

clean:
	rm -f $(objects) *~ *.out junk

cleanall: clean
	rm -f atlas

$(objects) : %.o : %.cpp
	$(cc) $(cflags) -o $*.o $*.cpp

$(dependencies) : %.d : %.cpp
	@$(cc) $(includedirs) -MM -MT $*.o $*.cpp

depend: $(dependencies)

# dependencies --- these were generated automatically by the command 
# make depend > make_dependencies on my system. Only local dependencies
# are considered. If you add new #include directives you should add the
# corresponding dependencies to the make_dependencies file; the best way if 
# your compiler supports the -MM option is probably to simply say again
# make depend > make_dependencies.

include make_dependencies

