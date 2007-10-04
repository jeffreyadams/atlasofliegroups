# See the file INSTALL for detailed instructions
# the following line avoids trouble on some systems (GNU make does this anyway)
SHELL = /bin/sh
INSTALL = /usr/bin/install

# You may edit the INSTALLDIR and BINDIR variables

# INSTALLDIR is where the executable atlas.exe and the messages directory
# will be moved after successfull compilation

# BINDIR is where a symbolic link 'atlas' to the executable will be placed;
# provided this is in your search path, you can then execute by typing: atlas

# In a single-user situation, you might want something like this:
#   INSTALLDIR := /home/fokko/myatlas
#   BINDIR     := /home/fokko/bin

# In a multi-user situation, you might want this (requires root):
#   INSTALLDIR := /usr/local/atlas
#   BINDIR     := /usr/local/bin

# The default is to use the current directory for both:
INSTALLDIR := $(shell pwd)
BINDIR := $(INSTALLDIR)

#Don't edit below this line, with the possible exception
#of rlincludes
###############################

MESSAGEDIR := $(INSTALLDIR)/messages/

# sourcedirs contains subdirectories of 'atlas/sources' that need compilation
sourcedirs := utilities error structure gkmod io interface test

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
 gflags := -c $(includedirs) -Wall -ggdb
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

# ensure that the software can find the help messages
cflags += -DMESSAGE_DIR_MACRO=\"$(MESSAGEDIR)\"

# the default compiler
CXX = g++

# give compiler=icc argument to make to use the intel compiler
ifdef compiler
    CXX = $(compiler)
endif

LDFLAGS := $(rlincludes)

# RULES follow below

# This target causes failed actions to clean up their (corrupted) target
.DELETE_ON_ERROR:

# we use no suffix rules
.SUFFIXES:

# The default target is 'all', which builds the executable
.PHONY: all
all: atlas.exe

# the following dependency forces emptymode.cpp to be recompiled whenever any
# of the object files changes; this guarantees that the date in the version
# string it prints will be that of the last recompilation.
sources/interface/emptymode.o: $(sources)

# For profiling not only 'cflags' used in compiling is modified, but linking
# also is different
atlas.exe: $(objects)
ifeq ($(profile),true)
	$(CXX) -pg -o atlas.exe $(objects) $(LDFLAGS)
else
	$(CXX) -o atlas.exe $(objects) $(LDFLAGS)
endif

install: atlas.exe
ifneq ($(INSTALLDIR),$(shell pwd))
	@echo "Installing directories and files in $(INSTALLDIR)"
	$(INSTALL) -d -m 755 $(INSTALLDIR)/www
	$(INSTALL) -d -m 755 $(INSTALLDIR)/messages
	$(INSTALL) -m 644 README $(INSTALLDIR)
	$(INSTALL) -m 755 atlas.exe $(INSTALLDIR)
	$(INSTALL) -m 644 www/*html $(INSTALLDIR)/www/
	$(INSTALL) -m 644 messages/*help $(INSTALLDIR)/messages/
	$(INSTALL) -m 644 messages/*intro_mess $(INSTALLDIR)/messages/
endif
	@if test -h $(BINDIR)/atlas; then rm -f $(BINDIR)/atlas; fi
	@if test -e $(BINDIR)/atlas ;\
	 then echo Warning: $(BINDIR)/atlas is not a symlink, I will not overwrite it;\
	 else ln -s $(INSTALLDIR)/atlas.exe $(BINDIR)/atlas ; fi

.PHONY: mostlyclean clean cleanall
mostlyclean:
	rm -f $(objects) *~ *.out junk

clean: mostlyclean
	rm -f atlas.exe

cleanall: clean
	rm -f $(dependencies)

realex:
	cd sources/interpreter && $(MAKE)

# The following two rules are static pattern rules: they are like implicit
# rules, but only apply to the files listed in $(objects) and $(dependencies).

$(objects) : %.o : %.cpp
	$(CXX) $(cflags) -o $*.o $*.cpp

$(dependencies) : %.d : %.cpp
	$(CXX) $(includedirs) -MM -MF $@ -MT "$*.o $*.d" $*.cpp

# now include all the files constructed by the previous rule
# this defines the dependencies, found by the preprocessor scanning #include
# directives, for all object files of the atlas
# make will automatically remake any of $(dependencies) if necessary

include $(dependencies)
