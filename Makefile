# a streamlined PMCAMx Makefile written by Jinhyok Heo <heoj@cmu.edu>, 2013
#
# run 'make C=case E=emission' to build

ifndef C
  C = PGH10
endif

ifndef E
  E = Base
endif

case = $(C)
emis = $(E)
app = PMCAMx_$(case)_$(emis)

# extract last two chacters as # of tags
tags = $(shell echo ${case}| awk '{print substr ($$0, length($$0)-1, 2)}')

srcExt = f
srcDir = src
buildDir = build
binDir = bin
inc = src/Inc

# compiler flags
#-O2             : set the optimization level to 2
#-tp k8-64       : target k8-64 AMD Opteron or Athlon-64 in 64-bit mode
#-pc 64          : set floating point precision to 64 double precision
#-Kieee          : in strict with the IEEE 754: less opt, more accurate math
#-Mdalign        : align doubles in structures on 8-byte boundaries
#-Mextend        : allow 132-column
#-Mnoframe       : don't set up a true stack frame pointer for functions
#-byteswapio     : swap bytes from small-endian to big-endian
#-mcmodel=medium : allows objects to be larger than 2GB => cannot be compiled statically

# Flags = -I$(inc) -O2 -tp k8-64 -pc 64 -Kieee -Mdalign -Mextend -Mnoframe -byteswapio -mcmodel=medium -L/usr/pgi/linux86-64/6.1/libso -lpgftnrtl
# Flags = -I$(inc) -O2 -tp k8-64 -pc 64 -Kieee -Mdalign -Mextend -Mnoframe -byteswapio -mcmodel=medium
# FC = /usr/pgi/linux86-64/6.1/bin/pgf77

# for Intel Fortran compiler
FC = ifort
Flags = -I$(inc) -fp-stack-check -g -traceback -heap-arrays -mcmodel=medium -shared-intel -O2 -align dcommons -extend_source -convert big_endian

sources := $(shell find $(srcDir) -name '*.$(srcExt)')
srcDirs := $(shell find $(srcDir) -type d -name .git -prune -o -type d)
objects := $(patsubst %.$(srcExt),$(buildDir)/%.o,$(sources))
# objects_July := $(patsubst %.$(srcExt),$(buildDir)/%.July.o,$(sources))


#all: $(binDir)/$(app) \
#     $(binDir)/$(app)_July
all: $(binDir)/$(app)
	@echo "Done!"

$(binDir)/$(app): buildrepo case $(objects)
	@echo "Linking $@"
	@$(FC) $(objects) $(Flags) -o $@

#$(binDir)/$(app)_July: buildrepo case_July emis_July $(objects_July)
#	@echo "Linking $@"
#	@$(FC) $(objects_July) $(Flags) -o $@

$(buildDir)/%.o: %.$(srcExt)
	@echo "Compiling $<"
	@$(FC) -c -o $@ $(Flags) $<

#$(buildDir)/%.July.o: %.$(srcExt)
#	@echo "Compiling $<"
#	@$(FC) -c -o $@ $(Flags) $<

case:
	ln -fs camx.prm.ladco.$(tags) $(srcDir)/Inc/camx.prm
	ln -fs Appemiss.f.$(case) $(srcDir)/APP/Appemiss.f
	ln -fs readar.f.$(emis) $(srcDir)/IO_bin/readar.f
	ln -fs readpt.f.$(emis) $(srcDir)/IO_bin/readpt.f

clean:
	$(RM) -r $(buildDir)

distclean: clean
	$(RM) -r $(binDir)/$(app)

buildrepo:
	@$(call make-repo)

# test:
# 	#$(FC) -o $(binDir)/test $(Flags) test.f
# 	echo $(tags)
# 	ln -fs camx.prm.ladco.$(tags) $(srcDir)/Inc/camx.prm

define make-repo
	for dir in $(srcDirs); do \
		mkdir -p $(buildDir)/$$dir; \
	done
endef
