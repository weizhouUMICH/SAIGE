
MAKEFILES_PATH := $(dir $(lastword $(MAKEFILE_LIST)))
include $(MAKEFILES_PATH)Makefile.base

RELEASE_FILE?=$(DIR_NAME).$(VERSION).tgz

ADDITIONAL_HELP= @echo "make install      Install binaries in $(INSTALLDIR)";\
	echo "make install INSTALLDIR=directory_for_binaries";\
	echo "                  Install binaries in directory_for_binaries"


.PHONY: package wholepackage

# Does not include the library.
package : 
# the touch gets rid of a tar warning
	touch $(RELEASE_FILE)
	tar chvz --exclude="*~" --exclude=$(RELEASE_FILE) --exclude='obj/*' --exclude='*.a'  --exclude='include/*' --exclude='bin/*' --exclude='test/results/*' --exclude-vcs -f $(RELEASE_FILE) --transform 's,^,$(DIR_NAME)_$(VERSION)/,' * --show-transformed-names 

BASE_LIB_PARTS := $(subst /, , $(BASE_LIB_PATH))
BASE_LIB_DIRNAME := $(word $(words $(BASE_LIB_PARTS)), $(BASE_LIB_PARTS))
WHOLEPACKAGE_MAKE := $(BASE_LIB_DIRNAME)/Makefiles/Makefile.wholepackage

DIR_ABOVE_LIB :=  $(patsubst %$(BASE_LIB_DIRNAME)/, %, $(BASE_LIB_PATH))

# also includes the library
wholepackage: 
# the touch gets rid of a tar warning
	touch $(RELEASE_FILE)
	tar chvz --exclude="*~" --exclude=$(RELEASE_FILE) --exclude='obj/*' --exclude='*.a'  --exclude='include/*' --exclude='bin/*' --exclude='test/results/*' --exclude-vcs -f $(RELEASE_FILE) --transform 's,^,$(DIR_NAME)_$(VERSION)/,;s,$(WHOLEPACKAGE_MAKE),Makefile,' -C .. $(DIR_NAME) -C $(DIR_NAME) -C $(DIR_ABOVE_LIB) $(BASE_LIB_DIRNAME) --show-transformed-names
