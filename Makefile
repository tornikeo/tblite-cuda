# CUDA Makefile

# Compiler and flags
NVCC := nvcc
CXXFLAGS := -O0 -std=c++20 -rdc=true --expt-relaxed-constexpr

# Directories
SRC_DIR := src
OBJ_DIR := obj
BIN_DIR := bin

# Files
TARGET := $(BIN_DIR)/main
SRCS := $(wildcard $(SRC_DIR)/*.cu)
OBJS := $(patsubst $(SRC_DIR)/%.cu, $(OBJ_DIR)/%.o, $(SRCS))

# Rules
all: $(TARGET)

$(TARGET): $(OBJS)
	mkdir -p $(BIN_DIR)
	$(NVCC) $(CXXFLAGS) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cu
	mkdir -p $(OBJ_DIR)
	$(NVCC) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

.PHONY: all clean