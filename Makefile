
# Compiler
CC = gcc

# Directories
SRC_DIR = src
INC_DIR = include
BUILD_DIR = build
FAST_DIR = build/fast

# Standard compiler flags
CFLAGS = -std=c11 -O2 -Wall -Wextra -Wpedantic
CPPFLAGS = -I$(INC_DIR) -D_POSIX_C_SOURCE=200809L
LDLIBS = -lm

# Fast compiler flags
FAST_CFLAGS = -std=c11 -O3 -march=native -DNDEBUG -Wall -Wextra -Wpedantic

# Sources
SRC = $(wildcard $(SRC_DIR)/*.c)

# Standard objects
OBJ = $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SRC))
DEP = $(OBJ:.o=.d)

# Fast objects
FAST_OBJ = $(patsubst $(SRC_DIR)/%.c,$(FAST_DIR)/%.o,$(SRC))
FAST_DEP = $(FAST_OBJ:.o=.d)

# Executable
TARGET = u1gauge2d

# Standard compilation
all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDLIBS)

# Fast compilation
fast: $(FAST_OBJ)
	$(CC) $(FAST_OBJ) -o $(TARGET) $(LDLIBS)

# Run
run: $(TARGET)
	./$(TARGET)

run-fast: fast
	./$(TARGET)

# Clean
clean:
	rm -rf $(BUILD_DIR)/*

.PHONY: all fast run run-fast clean
