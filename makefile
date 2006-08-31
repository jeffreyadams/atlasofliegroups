# the following line avoids trouble on some systems (GNU make does this anyway)
SHELL = /bin/sh
#where to put the wrapper script:
BINDIR = /usr/local/bin
#where the executable will be, change this line if this will be moved after compilation:
INSTALLDIR := $(shell pwd)
messagedir := $(INSTALLDIR)/messages/


# we use no suffix rules
.SUFFIXES:

# sourcedirs contains subdirectories of 'atlas/sources' that need compilation
sourcedirs := utilities memory error structure gkmod io interface test

# sources contains a list of the source files (i.e., the .cpp files)
sources := $(wildcard $(sourcedirs:%=sources/%/*.cpp))

# there is one .o file for each .cpp file, in the same directory
objects := $(sources:%.cpp=%.o)

# the variable 'dependencies' is used to automatically generate the
# dependencies, but the files described by this variable will not be made
dependencies := $(sources:%.cpp=%.d)

# headers are searched in the all directories containing source files
includedirs := $(addprefix -Isources/,$(sourcedirs))

# compiler flags
# whenever changing the setting of these flags, do 'make clean' first

# the following are predefined flavors for the compiler flags:
# optimizing (oflags), development with debugging (gflags), profiling (pflags)

oflags := -c $(includedirs) -Wall -O3 -DNDEBUG
gflags := -c $(includedirs) -Wall -g
pflags := -c $(includedirs) -Wall -pg -O -DNREADLINE

# the default setting is optimizing
cflags ?= $(oflags)

# to select another flavor, set debug=true or profile=true when calling make
# alternatively, you can set cflags="your personal flavor" to override default

ifeq ($(debug),true)
      cflags := $(gflags)
else
  ifeq ($(profile),true)
      cflags := $(pflags)
  endif
endif

# suppress readline by setting readline=false
# and/or make more verbose by setting verbose=true

ifeq ($(readline),false)
    cflags += -DNREADLINE
else
    rlincludes ?= -lreadline -lcurses
# to override this, either define and export a shell variable 'rlincludes'
# or set it when calling make. For instance for readline on the Mac do:
# $ make rlincludes="-lreadline.5.1 -lcurses"
endif

ifeq ($(verbose),true)
    cflags += -DVERBOSE
endif

CXX = g++ # the default compiler

# give compiler=icc argument to make to use the intel compiler
ifdef compiler
    CXX = $(compiler)
endif

# RULES follow below

# This target causes failed actions to clean up their (corrupted) target
.DELETE_ON_ERROR:

# The default target is 'all', which builds the executable and the wrapper
.PHONY: all
all: atlas

# For profiling not only 'cflags' used in compiling is modified, but linking
# also is different

cflags += -DMESSAGE_DIR_MACRO=\"$(messagedir)\"
atlas: $(objects)
ifeq ($(profile),true)
	$(CXX) -pg -o atlas.exe $(objects)
else
	$(CXX) -o atlas.exe $(objects) $(rlincludes)
endif
	./make_wrapper

install: 
	cp atlas $(BINDIR)

.PHONY: clean cleanall
clean:
	rm -f $(objects) *~ *.out junk

cleanall: clean
	rm -f atlas

# The following two rules are static pattern rules: they are like implicit
# rules, but only apply to the files listed in $(objects) and $(dependencies).
# However the files named in $(dependencies) are not actually made, and the
# call $(CXX) of g++ in its command only produces its dependency information
# on stdout (whence the prefix '@', to avoid the command itself on stdout)

$(objects) : %.o : %.cpp
	$(CXX) $(cflags) -o $*.o $*.cpp

$(dependencies) : %.d : %.cpp
	@$(CXX) $(includedirs) -MM -MT $*.o $*.cpp

.PHONY: depend
depend: $(dependencies)

# dependencies --- these were generated automatically by the command 
# make depend > make_dependencies on my system. Only local dependencies
# are considered. If you add new #include directives you should add the
# corresponding dependencies to the make_dependencies file; the best way if 
# your compiler supports the -MM option is probably to simply say again
# make depend > make_dependencies.

include make_dependencies



