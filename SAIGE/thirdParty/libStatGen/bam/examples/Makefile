PATH_TO_BASE=../../..
include ../../Makefiles/Makefile.include

SUBDIRS = $(wildcard */)

# Build in all subdirectories.
#
# see http://www.gnu.org/software/make/manual/make.html#Phony-Targets
# for a way of improving the following:
#

.PHONY : clean $(SUBDIRS)

all: $(SUBDIRS)

$(SUBDIRS):
	(echo "building in directory $@"; $(MAKE) $(PARALLEL_MAKE) OPTFLAG="$(OPTFLAG)" --no-print-directory -C $@ ${RECURSIVE_TARGET})

test: RECURSIVE_TARGET = test
test: $(SUBDIRS)

clean: RECURSIVE_TARGET = clean
clean: $(SUBDIRS)
