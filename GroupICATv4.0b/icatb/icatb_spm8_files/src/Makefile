#!/usr/bin/env make -f
# General Makefile to compile SPM C-MEX files
#
# Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
#
# $Id: Makefile 4199 2011-02-10 20:07:17Z guillaume $
#
###############################################################################
#
# This Makefile has been tested under Linux, Windows and MacOS.
# 
# If you have to tweak this Makefile or Makefile.var to compile the SPM 
# mex-files for your platform, please send the details to <spm@fil.ion.ucl.ac.uk> 
# so they can be included here. 
#
# To compile and install SPM, type the following in a Unix console:
# >  make && make install
#
# You can specify a particular platform with the following syntax:
# >  make PLATFORM=Your_Platform
# The standard targets are 'all', 'clean', 'distclean' and 'install'.
#
# For a list of compatible compilers, see
#    http://www.mathworks.com/support/compilers/current_release/
#
###############################################################################

include Makefile.var

###############################################################################
# Objects to go in the archive and mexfiles
###############################################################################

OBS     =\
	utils_uchar.$(SUF).o utils_short.$(SUF).o utils_int.$(SUF).o \
	utils_schar.$(SUF).o utils_ushort.$(SUF).o utils_uint.$(SUF).o\
	utils_float.$(SUF).o utils_double.$(SUF).o\
	utils_short_s.$(SUF).o utils_int_s.$(SUF).o\
	utils_ushort_s.$(SUF).o utils_uint_s.$(SUF).o\
	utils_float_s.$(SUF).o utils_double_s.$(SUF).o\
	spm_make_lookup.$(SUF).o spm_getdata.$(SUF).o spm_vol_access.$(SUF).o\
	spm_mapping.$(SUF).o

SPMMEX  =\
	icatb_spm_sample_vol.$(SUF) icatb_spm_slice_vol.$(SUF) icatb_spm_existfile.$(SUF) icatb_spm_bwlabel.$(SUF)\
	icatb_file2mat.$(SUF) icatb_mat2file.$(SUF)\

SUBDIRS =\
	toolbox/DARTEL \
	toolbox/FieldMap \
	toolbox/HDW \
	@xmltree/private

###############################################################################
# The main ways to run make
###############################################################################

all: verb.$(SUF) $(SPMMEX) verb.end

clean: verb.clean
	$(DEL) $(OBS)
ifeq (mex,$(SUF))
	$(DEL) $(subst .$(SUF),.o,$(SPMMEX)) spm_matfuns.o
endif

distclean: clean verb.distclean
	$(DEL) $(SPMMEX) spm_vol_utils.$(SUF).a

archive: spm_vol_utils.$(SUF).a

install: verb.install $(SPMMEX)
	$(COPY) $(SPMMEX) ..
	$(MOVE) ../icatb_file2mat.$(SUF) ../icatb_mat2file.$(SUF) ../@icatb_file_array/private

toolbox toolbox-all toolbox-clean toolbox-distclean toolbox-install toolbox-tarball:
	@@for p in $(SUBDIRS); do \
	echo "_____________________________________________________________";\
	echo "";\
	echo "          In " ../$$p; \
	echo "_____________________________________________________________";\
	echo ""	;\
	$(MAKE) -C ../$$p --no-print-directory $(subst -,,$(subst toolbox,,$@));\
	done

external external-all external-clean external-distclean external-install external-tarball: verb.external
	@ $(MAKE) -C ../external --no-print-directory $(subst -,,$(subst external,,$@))

tarball: all toolbox-all toolbox-tarball verb.archive
	$(TAR) cf spm_mex.tar $(SPMMEX)
	@@for p in $(SUBDIRS); do \
		$(TAR) Af spm_mex.tar  ../$$p/spm_mex.tar;\
	done
	$(ZIP) spm_mex.tar

###############################################################################
# Compile spm_vol_utils.c with various flags
###############################################################################

spm_vol_utils.$(SUF).a: $(OBS)
	$(DEL) $@
ifeq (win64,$(PLATFORM))
	$(AR)$@ $(OBS)
else
	$(AR) $@ $(OBS)
endif

UTILS = spm_vol_utils.c spm_make_lookup.h spm_getdata.h

utils_uchar.$(SUF).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_UNSIGNED_CHAR $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_short.$(SUF).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_SIGNED_SHORT $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_int.$(SUF).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_SIGNED_INT $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_schar.$(SUF).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_SIGNED_CHAR $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_ushort.$(SUF).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_UNSIGNED_SHORT $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_uint.$(SUF).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_UNSIGNED_INT $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_float.$(SUF).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_FLOAT $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_double.$(SUF).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_DOUBLE $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_short_s.$(SUF).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_SIGNED_SHORT -DSPM_BYTESWAP $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_int_s.$(SUF).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_SIGNED_INT -DSPM_BYTESWAP $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_ushort_s.$(SUF).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_UNSIGNED_SHORT -DSPM_BYTESWAP $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_uint_s.$(SUF).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_UNSIGNED_INT -DSPM_BYTESWAP $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_float_s.$(SUF).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_FLOAT -DSPM_BYTESWAP $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

utils_double_s.$(SUF).o: $(UTILS)
	$(MEX) -c spm_vol_utils.c -DSPM_DOUBLE -DSPM_BYTESWAP $(MEXEND)
	$(MOVE) spm_vol_utils.$(MOSUF) $@

###############################################################################
# Compile a few additional C routines for linking
###############################################################################

%.$(SUF).o : %.c
	$(MEX) -c $< $(MEXEND)
	$(MOVE) %.$(MOSUF) $@

spm_getdata.$(SUF).o: spm_getdata.c spm_getdata.h
	$(MEX) -c spm_getdata.c $(MEXEND)
	$(MOVE) spm_getdata.$(MOSUF) $@
	
spm_vol_access.$(SUF).o: spm_vol_access.c spm_vol_access.h spm_datatypes.h
	$(MEX) -c spm_vol_access.c $(MEXEND)
	$(MOVE) spm_vol_access.$(MOSUF) $@

spm_make_lookup.$(SUF).o: spm_make_lookup.c spm_make_lookup.h
	$(MEX) -c spm_make_lookup.c $(MEXEND)
	$(MOVE) spm_make_lookup.$(MOSUF) $@
	
spm_mapping.$(SUF).o: spm_mapping.c spm_mapping.h spm_vol_access.h spm_datatypes.h
	$(MEX) -c spm_mapping.c $(MEXEND)
	$(MOVE) spm_mapping.$(MOSUF) $@

###############################################################################
# Compile the mex files themselves
###############################################################################

%.$(SUF) : %.c
	$(MEX) $< $(MEXEND)

icatb_spm_sample_vol.$(SUF): spm_sample_vol.c spm_vol_utils.$(SUF).a\
		spm_mapping.h spm_vol_access.h
	$(MEX) spm_sample_vol.c spm_vol_utils.$(SUF).a -output icatb_spm_sample_vol.$(SUF) $(MEXEND)

icatb_spm_slice_vol.$(SUF): spm_slice_vol.c  spm_vol_utils.$(SUF).a\
		spm_mapping.h spm_vol_access.h
	$(MEX) spm_slice_vol.c  spm_vol_utils.$(SUF).a -output icatb_spm_slice_vol.$(SUF) $(MEXEND)

icatb_file2mat.$(SUF): file2mat.c
	$(MEX) file2mat.c -output icatb_file2mat.$(SUF) $(MEXEND)

icatb_mat2file.$(SUF): mat2file.c
	$(MEX) mat2file.c -output icatb_mat2file.$(SUF) $(MEXEND)

icatb_spm_bwlabel.$(SUF): spm_bwlabel.c
	$(MEX) spm_bwlabel.c -output icatb_spm_bwlabel.$(SUF) $(MEXEND)

icatb_spm_existfile.$(SUF): spm_existfile.c
	$(MEX) spm_existfile.c -output icatb_spm_existfile.$(SUF) $(MEXEND)


###############################################################################
# General messages
###############################################################################

verb.end:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        FINISHED"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.clean:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Deleting object (.o) files"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.distclean:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Deleting mex and archive (.a) files"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.install:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Installing mex files"
	@ echo "_____________________________________________________________"
	@ echo ""	

verb.archive:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Creating archive spm_mex.tar.gz"
	@ echo "_____________________________________________________________"
	@ echo ""	

###############################################################################
# Assorted architecture dependent messages
###############################################################################

verb.mexw32:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "         Windows compilation (32 bit)"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexw64:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Windows compilation (64 bit)"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexglx:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Linux compilation (x86-32)"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexa64:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Linux compilation (x86-64)"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexmac:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Unix compilation (MacOS X, PowerPC)"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexmaci:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Unix compilation (MacOS X, Intel 32 bit)"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexmaci64:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Unix compilation (MacOS X, Intel 64 bit)"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexsol:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Unix compilation (Solaris 32 bit)"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mexs64:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        Unix compilation (Solaris 64 bit)"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.mex:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        ${PLATFORM} compilation (`octave -v | head -n 1`)"
	@ echo "_____________________________________________________________"
	@ echo ""

verb.external:
	@ echo "_____________________________________________________________"
	@ echo ""
	@ echo "        In external"
	@ echo "_____________________________________________________________"
	@ echo ""
