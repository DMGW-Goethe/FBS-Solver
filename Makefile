CC := g++

AR := ar rv

# Comment to remove OMP support
OMP:=1

OBJ_DIR := build
INC_DIR := include
SRC_DIR := src

SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

PYTHON := python3

PYFBS_DIR := pyfbs
PYFBS_FILES := $(PYFBS_DIR)/cpyfbs.pxd $(PYFBS_DIR)/pyfbs.pyx

CPPFLAGS := -I$(INC_DIR)
CFLAGS   := -Wall -O3 -g -fPIC
LDFLAGS :=
LDLIBS :=


#DEBUG_PLOTTING:=1

ifdef DEBUG_PLOTTING
	CPPFLAGS:=$(CPPFLAGS) -I/usr/include/python3.8 -DDEBUG_PLOTTING
	LDLIBS:=$(LDLIBS) -lpython3.8
endif

ifdef OMP
	CPPFLAGS:=$(CPPFLAGS) -fopenmp
endif

.PHONY: all clean archive

all: libfbs.a fbs pyfbs

fbs: $(OBJ)
	$(CC) $(CFLAGS) $(CPPFLAGS) $^ main.cpp -o main.out $(LDFLAGS) $(LDLIBS)

archive: libfbs.a

pyfbs: archive $(PYFBS_FILES)
	cd $(PYFBS_DIR); $(PYTHON) setup.py build_ext --inplace

libfbs.a: $(OBJ)
	$(AR)  $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.hpp | $(OBJ_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(OBJ_DIR):
	mkdir -p $@

clean:
	@$(RM) -rv  $(OBJ_DIR)
	@$(RM) main.out
	@$(RM) libfbs.a
	@$(RM) -v $(PYFBS_DIR)/pyfbs.cp*
