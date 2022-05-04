INSTALL?=cp
DESTDIR?=/usr
bindir?=/bin

install: $(programs:%=install_%)

install_%: $(DRAGEN_OS_BUILD)/%
	$(INSTALL) $<  $(DESTDIR)$(bindir)/$*
	
