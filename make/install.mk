INSTALL?=cp
DESTDIR?=/usr
bindir?=/bin

install: $(programs:%=install_%)

install_%: $(BUILD)/%
	$(INSTALL) $<  $(DESTDIR)$(bindir)/$*
	
