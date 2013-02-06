#
# Paul H J Kelly Imperial College London 2013
#
# Makefile structure from:
#  http://www.gnu.org/software/make/manual/html_node/make_44.html#SEC48
#
# nucleic.h is copied to /tmp to avoid a problem that the simplescalar
# compiler's preprocessor has a 32-bit/64-bit compatibility problem
# and can't handle includes of files on filesystems with 64-bit stats.


SRC=nucleic-core.c nucleic-database.c init_nucleotides.c

all: 	nucleic.x86 nucleic.ss


################################################################
# Native linux x86

nucleic.x86: $(SRC)
	cp -f nucleic.h /tmp/nucleic.h
	gcc -o $@ -O3 -g $(MYFLAGS) $(SRC) -lm


################################################################
# Simplescalar

SSGCC=/homes/phjk/simplescalar/bin/gcc

OBJss=$(SRC:.c=.ss)

simplescalar: $(OBJss) Makefile

nucleic.ss: $(SRC)
	cp -f nucleic.h /tmp/nucleic.h
	$(SSGCC) -o $@ -O3 -g $(MYFLAGS) $(SRC)


################################################################
# This rule has to come last; it deletes all the binaries.
# If you define a new class of binaries, add it to this list
#
clean:	
	rm -f nucleic.x86 nucleic.ss

.PHONY: all clean

