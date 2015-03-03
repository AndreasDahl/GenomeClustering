# Compiler flags
CXX			=g++
CPPFLAGS	+=-std=c++11 #-Wall
#CPPFLAGS	+=-static-libgcc -static-libstdc++ 
CPPFLAGS	+=-ansi -pedantic 

# For Debug
#CPPFLAGS	+=-Wextra
CPPFLAGS	+=-g

# For Optimization
CPPFLAGS	+=-O4

#ifeq ($(OS),Windows_NT)
#else
#endif

# File objects
SOURCES		=$(wildcard *.cpp) $(wildcard */*.cpp)
OBJECTS		=$(SOURCES:.cpp=.o)
WINOBJECTS	=$(subst /,\,$(OBJECTS))
EXECUTABLE	=dmsearch

.PHONY: default clean cleanexe

default: cleanexe $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
		$(CXX) $(CPPFLAGS) -o  $(EXECUTABLE) $(OBJECTS) $(LIBS) $(LDFLAGS)

clean: cleanexe
ifeq ($(OS),Windows_NT)
	del $(WINOBJECTS)
else
	rm -rf $(OBJECTS)
endif

cleanexe:
ifeq ($(OS),Windows_NT)
	del $(EXECUTABLE)
else
	rm -rf $(EXECUTABLE)
endif
