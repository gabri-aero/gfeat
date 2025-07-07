# Define compiler and flags
CXX = g++
CXX_FLAGS = -O3 -ffast-math -march=native -Wall -fPIC
LD_FLAGS = -shared

# Add data dir to flags
DATA_DIR = data
CXX_FLAGS += -DDATA_DIR=\"$(DATA_DIR)\"

# Define directories
INCLUDE_DIR = gfeat
BUILD_DIR = build
TEST_DIR = tests

# External directories
FFT_DIR = external/fft
FUNCTIONS_DIR = external/functions

# External libraries
EXTERNAL_INCLUDES = -I$(FFT_DIR) -I$(FUNCTIONS_DIR)

# Define header files and object files
HEADERS = $(shell find $(INCLUDE_DIR) -name '*.hpp')
TEST_SOURCES = $(shell find $(TEST_DIR) -name '*.cpp')
TEST_EXES = $(patsubst $(TEST_DIR)/%.cpp,$(BUILD_DIR)/%.exe,$(TEST_SOURCES))

# Define GTest libraries
GTEST_LIBS = -lgtest -lgtest_main -pthread
GTEST_DIR = /usr/include/gtest

# Build tests
test: $(TEST_EXES)
	@for test_exe in $(TEST_EXES); do \
		./$$test_exe; \
	done

# Create build directory if it does not exist
$(BUILD_DIR):
	@mkdir -p $(BUILD_DIR)

# Test targets
$(BUILD_DIR)/%.exe: $(TEST_DIR)/%.cpp $(BUILD_DIR) $(DATA_DIR)
	mkdir -p $(shell dirname $@)
	$(CXX) $(CXX_FLAGS) $(EXTERNAL_INCLUDES) -I. -I$(GTEST_DIR) $< -o $@ $(GTEST_LIBS)

# Download gravity field data
$(DATA_DIR):
	bash default_data.sh $@

# Clean build files
clean:
	rm -rf $(BUILD_DIR)
