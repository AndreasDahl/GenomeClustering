# Compiler flags
CXX			=g++
CPPFLAGS	+=-std=c++11
CPPFLAGS	+=-Wall
#CPPFLAGS	+=-static-libgcc -static-libstdc++

# For Debug
CPPFLAGS	+=-Wextra -Werror
#CPPFLAGS	+=-g

# For Optimization
CPPFLAGS	+=-O4

#ifeq ($(OS),Windows_NT)
#else
#endif

# File objects
SOURCES		=$(wildcard *.cpp) $(wildcard */*.cpp)
OBJECTS		=$(SOURCES:.cpp=.o)
WINOBJECTS	=$(subst /,\,$(OBJECTS))
OUTFOLDER	=out
EXECUTABLE	=$(OUTFOLDER)/dmsearch
WINEXE		=$(subst /,\,$(EXECUTABLE)).exe

.PHONY: default clean cleanexe createfolder

default: createfolder clean $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
		$(CXX) $(CPPFLAGS) -o  $(EXECUTABLE) $(OBJECTS) $(LIBS) $(LDFLAGS)

createfolder:
ifeq ($(OS),Windows_NT)
	if not exist "$(OUTFOLDER)" mkdir $(OUTFOLDER)
else
	mkdir -p $(OUTFOLDER)
endif

clean: cleanexe
ifeq ($(OS),Windows_NT)
	del $(WINOBJECTS)
else
	rm -rf $(OBJECTS)
endif

cleanexe:
ifeq ($(OS),Windows_NT)
	del $(WINEXE)
else
	rm -rf $(EXECUTABLE)
endif
