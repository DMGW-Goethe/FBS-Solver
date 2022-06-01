CC := g++

# Comment to remove OMP support
OMP:=1

OBJ_DIR :=build
INC_DIR :=include
SRC_DIR :=src

SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

CPPFLAGS := -I$(INC_DIR)
CFLAGS   := -Wall -O3 -g
LDFLAGS :=
LDLIBS :=


ifdef DEBUG_PLOTTING
	CPPFLAGS:=$(CPPFLAGS) -I/usr/include/python3.8 -DDEBUG_PLOTTING
	LDLIBS:=$(LDLIBS) -lpython3.8
endif

ifdef OMP
	CPPFLAGS:=$(CPPFLAGS) -fopenmp
endif

.PHONY: all clean

all: $(OBJ)
	echo $(OBJ)
	$(CC) $(CFLAGS) $(CPPFLAGS) $^ main.cpp -o main.out $(LDFLAGS) $(LDLIBS)


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.hpp | $(OBJ_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(OBJ_DIR):
	mkdir -p $@

clean:
	@$(RM) -rv  $(OBJ_DIR) main.out
