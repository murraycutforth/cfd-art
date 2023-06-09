#
# Makefile for the 2D-Allaire-diffuse project
#
# External dependencies
#	- The header-only Eigen library. Download from eigen.tuxfamily.org
#	- My own PhD-2D-HCLframework headers. Download from https://github.com/murraycutforth/PhD-2D-HCLsolver-framework
#	- My own PhD-Common headers. Download from https://github.com/murraycutforth/PhD-Common
#
# Based on the file at http://hiltmon.com/
#

# Target
TARGET := cfd-art

# Folders
BUILDDIR := objectfiles
SRCDIR := sourcecode
INCDIRS := headers
INCLIST := -I./headers

# Files
SOURCES := $(shell find $(SRCDIR) -type f -name *.cpp)
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.cpp=.o))

# Folder Lists
BUILDLIST := $(BUILDDIR)
EXTINCLIST := $(HCLFRAMEWORK) $(COMMON) $(INCLIST)

# Shared Compiler Flags
OPLEVEL := -O3
CFLAGS := -Wall -c -fopenmp -std=c++11 $(OPLEVEL)
LINKFLAGS := -fopenmp $(OPLEVEL)

# Linking Step
$(TARGET): $(OBJECTS)
	@mkdir -p $(BUILDLIST)
	@echo "Linking..."
	@echo "Linking $(TARGET) using options: $(LINKFLAGS)"; g++ $^ $(LINKFLAGS) -o $(TARGET)
	@echo "Success!"

# Compilation Step
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(BUILDLIST)
	@echo "Compiling $< using options: $(CFLAGS)"; g++ $(CFLAGS) $(EXTINCLIST) -o $@ $<

clean:
	@echo "Cleaning $(TARGET)..."; rm $(BUILDDIR)/* $(TARGET)
