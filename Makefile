TOPDIR=/Users/jinmei/software/build/libint
CWD=$(shell pwd)
ifndef SRCDIR
  SRCDIR=/Users/jinmei/software/source/libint
endif
-include $(TOPDIR)/src/bin/MakeVars
-include $(TOPDIR)/src/lib/libint/MakeVars.features

# include headers the object include directory
CPPFLAGS += -I$(TOPDIR)/include -I$(TOPDIR)/include/libint2 -I$(SRCDIR)/$(TOPDIR)/src/lib/libint -I/Users/jinmei/software/install/libint/include -I/usr/local/include/eigen3 -DSRCDATADIR=\"$(SRCDIR)/lib/basis\"

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


$(INTS): $(CXXINTSOBJ) $(COMPILER_LIB) $(COMPUTE_LIB)
	$(LD) -o $@ $(LDFLAGS) $^ $(SYSLIBS) 

# Source files for timer and tester are to be compiled using CXXGEN
$(INTS) : CXX=$(CXXGEN)
$(INTS) : CXXFLAGS=$(CXXGENFLAGS)
$(INTS) : LD=$(CXXGEN)

clean::
	-rm -rf $(INTS) *.o *.d

distclean:: realclean
	-rm -rf $(TOPDIR)/include/libint2/boost

realclean:: clean

targetclean:: clean

$(TOPDIR)/include/libint2/boost/preprocessor.hpp: $(SRCDIR)/$(TOPDIR)/external/boost.tar.gz
	gunzip -c $(SRCDIR)/$(TOPDIR)/external/boost.tar.gz | tar -xf - -C $(TOPDIR)/include/libint2

depend:: $(CXXINTSDEP) $(CXXTEST1DEP) $(CXXTEST2DEP)

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

