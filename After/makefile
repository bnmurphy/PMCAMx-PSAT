#  Makefile for After post-processing programs
#
#  Modified by Ben Murphy - 3/11/08
#
#################################################################
#
#  Portland Group Compiler Flags (pgf77)
#
#    -tp p6		specify the target processor 
#			  (p6=pentium pro/pentium II)
#    -pc 64		precision option (64=double precision)
#    -Kieee		Perform floating-point operations 
#			  in strict conformance with the IEEE 
#			  754 standard
#    -Mdalign		Align doubles in structures on 8-byte
#			  boundaries
#    -Mextend		Allow 132 column source lines
#    -Mnoframe		Do not set up a true stack frame pointer 
#			  for fxns
#    -byteswapio	Swap from big-endian to little-endian
#			  or vice-versa on input/output of
#			  unformatted FORTRAN data
#    -W1,		Pass each comma-delimeted option to 
#			  the linker
#       -Bstatic	Static linking
#    -Minfo		Emit useful info to stderr
#
#    -c			skip the link step; compile and assemble 
#			  only
#    -o			output file
#
#  Debugging Options
#
#    -g			generate symbolic debug info.
#			sets optimizations level to zero.
#
################################################################


FC = /usr/pgi/linux86-64/6.0/bin/pgf77

FLGS =  -pc 64 -Kieee -Mdalign -mcmodel medium -Mextend -Mnoframe -byteswapio -Wl,-Bstatic -Minfo -Minform=warn
#FLGS = -fast

AMB = cmp_improve.f cmp_airs.f cmp_paqs.f sub_after.f sub_after_di.f
AFT = after_saprc.f after_di_saprc.f after_di_saprc_vbins.f
AFTL = after_saprc_lay.f after_di_saprc.f after_di_saprc_vbins.f

App_After_Ben: App_After_Ben.f
	$(FC) $(FLGS) -c -g App_After_Ben.f
	$(FC) $(FLGS) -g App_After_Ben.o -o App_After_Ben


clean:
	rm *.o avgrads
