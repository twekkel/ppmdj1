## Take a look at PPMdType.h for additional compiler & environment options
PROJECT = PPMd
DEBUG = 0
CPP_SET = $(PROJECT).cpp Model.cpp
C_SET =

CC = gcc
LINK = gcc
CODE_FLAGS = -march=native -fno-exceptions -fno-rtti -pedantic -Wall \
	-Wno-unknown-pragmas
LIBS = 	

ifeq ($(DEBUG),0)
  OPT_FLAGS = -O3 -ffast-math -funroll-all-loops -floop-optimize
  DEBUG_FLAGS = -g0
  LD_FLAGS = -s -O3
else
  OPT_FLAGS = -O0
  DEBUG_FLAGS = -g
  LD_FLAGS = -g
endif

OBJ_SET = $(CPP_SET:.cpp=.o) $(C_SET:.c=.o)

$(PROJECT): $(OBJ_SET)
	@$(LINK) -o $(PROJECT) $(OBJ_SET) $(LIBS) -lstdc++ -Xlinker $(LD_FLAGS)

.cpp.o:
	@$(CC) $(CODE_FLAGS) $(OPT_FLAGS) $(DEBUG_FLAGS) -c $^
.c.o:
	@$(CC) $(CODE_FLAGS) $(OPT_FLAGS) $(DEBUG_FLAGS) -c $^
