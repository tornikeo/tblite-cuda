# Thanks to Job Vranish (https://spin.atomicobject.com/2016/08/26/makefile-c-projects/)
TARGET_EXEC := main.bin
TEST_EXEC := test.bin
NVCC := nvcc -ccbin=gcc -G -g -rdc=true -std=c++17
BUILD_DIR := ./build
SRC_DIRS := ./src

# Find all the C and C++ files we want to compile
# Note the single quotes around the * expressions. The shell will incorrectly expand these otherwise, but we want to send the * directly to the find command.
SRCS := $(shell find $(SRC_DIRS) -name '*.cu')

TEST_SRCS := $(shell find ./test -name '*.cu')
TEST_BINS := $(patsubst ./test/%.cu, $(BUILD_DIR)/test/%.bin, $(TEST_SRCS))

# Prepends BUILD_DIR and appends .o to every src file
# As an example, ./your_dir/hello.cpp turns into ./build/./your_dir/hello.cpp.o
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

# String substitution (suffix version without %).
# As an example, ./build/hello.cpp.o turns into ./build/hello.cpp.d
DEPS := $(OBJS:.o=.d)

# Every folder in ./src will need to be passed to GCC so that it can find header files
INC_DIRS := $(shell find $(SRC_DIRS) -type d)
# Add a prefix to INC_DIRS. So moduleA would become -ImoduleA. GCC understands this -I flag
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

# The -MMD and -MP flags together generate Makefiles for us!
# These files will have .d instead of .o as the output.
CPPFLAGS := -G -g $(INC_FLAGS)  -MMD -MP -DCUDA_SEPARABLE_COMPILATION
# Add the CUDA library path and link the cusolver library
CUDA_LIB := /home/tornikeo/micromamba/envs/pb/lib
export LD_LIBRARY_PATH := $(CUDA_LIB):$(LD_LIBRARY_PATH)
LDFLAGS := -L$(CUDA_LIB) -lcusolver -lcublas
NVCCFLAGS := -w -Xcompiler -Wno-unused-variable

# Run all test executables
test: $(TEST_BINS)
	for exe in $(TEST_BINS); do \
		echo "Running $$exe..."; \
		./$$exe || exit 1; \
	done

# Build each test executable from its corresponding test CUDA source file
$(BUILD_DIR)/test/%.bin: ./test/%.cu
	mkdir -p $(dir $@)
	$(NVCC) $< -o $@ $(CPPFLAGS) $(NVCCFLAGS)  $(LDFLAGS)

# The final build step.
$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	$(NVCC) $(OBJS) -o $@ $(CPPFLAGS) $(NVCCFLAGS) $(LDFLAGS)

# Build step for CUDA source
$(BUILD_DIR)/%.cu.o: %.cu
	mkdir -p $(dir $@)
	$(NVCC) -c $< -o $@ $(CPPFLAGS) $(NVCCFLAGS)

.PHONY: clean
clean:
	rm -r $(BUILD_DIR)

# Include the .d makefiles. The - at the front suppresses the errors of missing
# Makefiles. Initially, all the .d files will be missing, and we don't want those
# errors to show up.
-include $(DEPS)