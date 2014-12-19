#
# Makefile for a SVL test program
#

NAME        =   svltest

SRCS        =   SVLTest.cpp

LIBS        =   -lsvl -lm
DBG_LIBS    =   -lsvl.dbg -lm

INC_DIRS    =   -I../include
LIB_DIRS    =   -L../lib

BUILD_DBG_FLAGS = -DVL_DEBUG

include ../makefiles/prog.mf

test:	svltest.dbg
	-./svltest.dbg | \
	  sed -e 's/[0-9].[0-9]*e-1[6-9]/0/g; s/-0\([^\.0-9]\)/0\1/g' > out
	@-echo \
	  "--- Differences from reference run ------------------------------"
	@-diff out out.ref

# the sed command above clips numbers like 6.12345e-16 to zero, and converts
# -0 to 0. 
