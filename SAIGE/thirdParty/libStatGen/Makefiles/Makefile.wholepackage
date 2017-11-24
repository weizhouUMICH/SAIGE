SUBDIRS = $(sort $(dir $(wildcard */)))

SUBDIRS_NO_STATGEN = $(filter-out libStatGen/,$(SUBDIRS))

MAKEFILES_PATH := libStatGen/Makefiles/
include $(MAKEFILES_PATH)Makefile.include

# Build in all subdirectories.


.PHONY: $(SUBDIRS)

include $(MAKEFILES_PATH)Makefile.help

$(SUBDIRS): 
	@$(MAKE) -C $@ $(MAKECMDGOALS)

$(SUBDIRS_NO_STATGEN): libStatGen/

%: $(SUBDIRS) ;

Makefile.%: ;

Makefile: ;
SUBDIRS = $(sort $(dir $(wildcard */)))
