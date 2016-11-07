TOPDIR= libint build directory (/Users/jinmei/software/build/libint)
INSDIR= libint install directory (/Users/jinmei/software/install/libint)
CWD=$(shell pwd)
ifndef SRCDIR
  SRCDIR=libint source directory (/Users/jinmei/software/source/libint)
endif
-include $(TOPDIR)/src/bin/MakeVars
-include $(TOPDIR)/src/lib/libint/MakeVars.features

# include headers the object include directory
CPPFLAGS += -I$(TOPDIR)/include -I$(TOPDIR)/include/libint2 -I$(SRCDIR)/$(TOPDIR)/src/lib/libint -I/($INSDIR)/include -I/usr/local/include/eigen3 -DSRCDATADIR=\"$(SRCDIR)/lib/basis\"

COMPILER_LIB = $(TOPDIR)/src/bin/libint/libINT.a
COMPUTE_LIB = -lint2
vpath %.a $(TOPDIR)/lib:$(TOPDIR)/lib/.libs

OBJSUF = o
DEPSUF = d
CXXDEPENDSUF = none
CXXDEPENDFLAGS = -M

INTS = ints_2e
CXXINTSSRC = $(INTS).cc
CXXINTSOBJ = $(CXXINTSSRC:%.cc=%.$(OBJSUF))
CXXINTSDEP = $(CXXINTSSRC:%.cc=%.$(DEPSUF))

TEST1 = hartree-fock
CXXTEST1SRC = $(TEST1).cc
CXXTEST1OBJ = $(CXXTEST1SRC:%.cc=%.$(OBJSUF))
CXXTEST1DEP = $(CXXTEST1SRC:%.cc=%.$(DEPSUF))

TEST2 = hartree-fock++
CXXTEST2SRC = $(TEST2).cc
CXXTEST2OBJ = $(CXXTEST2SRC:%.cc=%.$(OBJSUF))
CXXTEST2DEP = $(CXXTEST2SRC:%.cc=%.$(DEPSUF))


$(INTS): $(CXXINTSOBJ) $(COMPILER_LIB) $(COMPUTE_LIB)
	$(LD) -o $@ $(LDFLAGS) $^ $(SYSLIBS) 
	
$(TEST1): $(CXXTEST1OBJ) $(COMPILER_LIB) $(COMPUTE_LIB)
	$(LD) -o $@ $(LDFLAGS) $^ $(SYSLIBS)

$(TEST2): $(CXXTEST2OBJ) $(COMPILER_LIB) $(COMPUTE_LIB)
	$(LD) -o $@ $(LDFLAGS) $^ $(SYSLIBS)

# Source files for timer and tester are to be compiled using CXXGEN
$(INTS) $(TEST1) $(TEST2): CXX=$(CXXGEN)
$(INTS) $(TEST1) $(TEST2): CXXFLAGS=$(CXXGENFLAGS)
$(INTS) $(TEST1) $(TEST2): LD=$(CXXGEN)

clean::
	-rm -rf $(INTS) $(TEST1) $(TEST2) *.o *.d

distclean:: realclean
	-rm -rf $(TOPDIR)/include/libint2/boost

realclean:: clean

targetclean:: clean

$(TOPDIR)/include/libint2/boost/preprocessor.hpp: $(SRCDIR)/$(TOPDIR)/external/boost.tar.gz
	gunzip -c $(SRCDIR)/$(TOPDIR)/external/boost.tar.gz | tar -xf - -C $(TOPDIR)/include/libint2

depend:: $(CXXINTSDEP) $(CXXTEST1DEP)  $(CXXTEST2DEP) 

ifneq ($(CXXDEPENDSUF),none)
%.d:: %.cc $(TOPDIR)/include/libint2/boost/preprocessor.hpp
	$(CXXDEPEND) $(CXXDEPENDFLAGS) -c $(CPPFLAGS) $(CXXFLAGS) $< > /dev/null
	sed 's/^$*.o/$*.$(OBJSUF) $*.d/g' < $(*F).$(CXXDEPENDSUF) > $(@F)
	/bin/rm -f $(*F).$(CXXDEPENDSUF)
else
%.d:: %.cc $(TOPDIR)/include/libint2/boost/preprocessor.hpp
	$(CXXDEPEND) $(CXXDEPENDFLAGS) -c $(CPPFLAGS) $(CXXFLAGS) $< | sed 's/^$*.o/$*.$(OBJSUF) $*.d/g' > $(@F)
endif

-include $(CXXINTSDEP)
-include $(CXXTEST1DEP)
-include $(CXXTEST2DEP) 
