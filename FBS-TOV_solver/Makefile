CC := g++


OBJ_DIR :=build
INC_DIR :=include
SRC_DIR :=src

SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

CPPFLAGS := -I/usr/include/python3.8  -I$(INC_DIR)
CFLAGS   := -Wall -O3 -g
LDFLAGS :=
LDLIBS = -lpython3.8


.PHONY: all clean

all: $(OBJ)
	echo $(OBJ)
	$(CC) $(CFLAGS) $(CPPFLAGS) $^ main.cpp -o main $(LDFLAGS) $(LDLIBS)


$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(INC_DIR)/%.hpp | $(OBJ_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(OBJ_DIR):
	mkdir -p $@

clean:
	@$(RM) -rv  $(OBJ_DIR) main
