# project name 
TARGET	= trimFilter


#---------------------------------------------------------------------
# executables 
#---------------------------------------------------------------------
C	:= gcc 
MD	:= mkdir 
LINKER  := gcc -o
RM      := rm -f

# compiling flags
CFLAGS   = -Wall -O3 -march=native -std=c11  -Iinclude/  
LFLAGS   = -Wall   



# change these to set the proper directories where each files should be
SRCDIR   = src
INCDIR   = include
OBJDIR   = obj
BINDIR   = bin

SOURCES  := $(wildcard $(SRCDIR)/*.c)
INCLUDES := $(wildcard $(INCDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.c=$(OBJDIR)/%.o)


#---------------------------------------------------------------------
# rules
#---------------------------------------------------------------------

all: mkdirs $(BINDIR)/$(TARGET) 

mkdirs: 
	$(MD) -p $(OBJDIR)
	$(MD) -p $(BINDIR)

$(BINDIR)/$(TARGET): $(OBJECTS)
	@ echo "Linking..."
	$(LINKER) $@  $(LFLAGS)   $(OBJECTS)
	@echo "Linking complete!"

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.c
	@echo "Compiling ..."
	$(C) $(CFLAGS)  -c $< -o $@
	@echo "Compiled "$<" successfully!"


.PHONY: clean
clean:
	$(RM) $(OBJECTS)
	@echo "Cleanup complete!"

.PHONY: remove
remove: clean
	$(RM) $(BINDIR)/$(TARGET)
	@echo "Executable removed!"
 
