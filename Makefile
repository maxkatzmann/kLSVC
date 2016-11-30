CC := g++ # This is the main compiler
SRCDIR := src
BUILDDIR := build
TARGET := bin/kLSVC
CFLAGS := -std=c++11 -O3 $(DEBUG) $(REDUCEOUTPUT) -w -Wall
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name "*.$(SRCEXT)")
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

INC := -I include

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

debug:
	$(MAKE) $(MAKEFILE) DEBUG="-DDEBUG"

reducedoutput:
	$(MAKE) $(MAKEFILE) REDUCEOUTPUT="-DREDUCEDOUTPUT"

.PHONY: clean
