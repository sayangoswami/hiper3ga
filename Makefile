
# By default, it builds for ibv conduit. If you use a different conduit, export CONDUIT=<name> before running make
CONDUIT?=ibv
$(info Building for $(CONDUIT) conduit)

# Replace the GASNet installation path with the appropriate one
MAK=/usr/local/gasnet/include/$(CONDUIT)-conduit/$(CONDUIT)-par.mak
$(info Including file $(MAK))
include $(MAK)

CC         := $(GASNET_CC)
CXX        := $(GASNET_CXX)
LD         := $(GASNET_LD)

TARGET     := hiper3ga
SRCDIR     := src
INCDIR     := include
BUILDDIR   := obj
TARGETDIR  := debug
SRCEXT     := c
DEPEXT     := d
OBJEXT     := o

# Replace the libvbyte installation path with the appropriate one
CFLAGS     := $(GASNET_CFLAGS) -fdiagnostics-color=always
LDFLAGS    := $(GASNET_LDFLAGS) -L/usr/local/libvbyte
LIB        := $(GASNET_LIBS) -lvbyte -lm
INC        := $(GASNET_CPPFLAGS) -I/usr/local/libvbyte -I$(INCDIR) -DGASNETT_THREAD_SAFE=1
INCDEP     := -I$(INCDIR)

#---------------------------------------------------------------------------------
#DO NOT EDIT BELOW THIS LINE
#---------------------------------------------------------------------------------

SOURCES    := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS    := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))

#Defauilt Make
all: resources $(TARGET)

#Remake
remake: cleaner all

#Copy Resources from Resources Directory to Target Directory
resources: directories

#Make the Directories
directories:
	@mkdir -p $(TARGETDIR)
	@mkdir -p $(BUILDDIR)

#Clean only Objecst
clean:
	@$(RM) -rf $(BUILDDIR)

#Full Clean, Objects and Binaries
cleaner: clean
	@$(RM) -rf $(TARGETDIR)

#Pull in dependency info for *existing* .o files
-include $(OBJECTS:.$(OBJEXT)=.$(DEPEXT))

#Link executable
$(TARGET): $(OBJECTS)
	$(LD) $(LDFLAGS) -o $(TARGETDIR)/$(TARGET) $^ $(LIB)

#Compile
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<
	@$(CC) $(CFLAGS) $(INCDEP) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
	@rm -f $(BUILDDIR)/$*.$(DEPEXT).tmp