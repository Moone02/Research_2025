## Top
## Makefile
# Compiler and Linker
CC = gcc
# CFLAGS: Compiler flags
# -Wall -Wextra: Enable many useful warnings
# -std=c11: Use the C11 standard
# -g: Include debugging information
# -O2: Optimization level (optional, can remove for easier debugging initially)
CFLAGS = -Wall -Wextra -std=c11 -g -O2 
# LDFLAGS: Linker flags
# -lm: Link math library
# -lgsl -lgslcblas: Link GSL and GSL's CBLAS implementation
LDFLAGS = -lm -lgsl -lgslcblas

# --- Adjust GSL paths if not in standard locations ---
# Example for Homebrew on ARM Mac (Apple Silicon):
# GSL_INCLUDE_DIR = /opt/homebrew/opt/gsl/include
# GSL_LIB_DIR = /opt/homebrew/opt/gsl/lib
# Example for Homebrew on Intel Mac:
# GSL_INCLUDE_DIR = /usr/local/opt/gsl/include
# GSL_LIB_DIR = /usr/local/opt/gsl/lib
# If GSL_INCLUDE_DIR and GSL_LIB_DIR are set, uncomment the following:
# CFLAGS += -I$(GSL_INCLUDE_DIR)
# LDFLAGS += -L$(GSL_LIB_DIR)
# --- End GSL Path Adjustment ---

# Header file (dependency for most .c files)
HEADER = schwarzschild_tracer.h

# Source files for the library part (excluding main_test.c)
SRCS_LIB = \
    trajectory_module.c \
    image_mapping_utils.c \
    results_generation.c \
    photon_rendering.c \
    photon_chunk_io.c \
    ppm_image_io.c

# Object files for the library part
OBJS_LIB = $(SRCS_LIB:.c=.o)

# Main test program
MAIN_TEST_EXEC = schwarzschild_test
MAIN_TEST_SRC = main_test.c
MAIN_TEST_OBJ = $(MAIN_TEST_SRC:.c=.o)

# Default target: build the main test executable
all: $(MAIN_TEST_EXEC)

# Rule to link the main test executable
$(MAIN_TEST_EXEC): $(OBJS_LIB) $(MAIN_TEST_OBJ)
	@echo "Linking $(MAIN_TEST_EXEC)..."
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)
	@echo "$(MAIN_TEST_EXEC) built successfully."

# Generic rule to compile .c files into .o object files
# Depends on the source file itself and the main header
%.o: %.c $(HEADER)
	@echo "Compiling $<..."
	$(CC) $(CFLAGS) -c $< -o $@

# --- Unit Test Targets ---
# For each .c file that has an embedded main() for unit testing,
# create a target. The executable will be named "test_MODULENAME".
# The -DUNIT_TEST_MODULENAME flag activates the embedded main().

# Unit test for trajectory_module.c
UT_TRAJ_EXEC = test_trajectory_module
UT_TRAJ_SRC = trajectory_module.c
$(UT_TRAJ_EXEC): $(UT_TRAJ_SRC) $(HEADER)
	@echo "Building unit test $(UT_TRAJ_EXEC)..."
	$(CC) $(CFLAGS) -DUNIT_TEST_TRAJECTORY_MODULE $(UT_TRAJ_SRC) -o $@ $(LDFLAGS)
	@echo "Running unit test $(UT_TRAJ_EXEC)..."
	./$(UT_TRAJ_EXEC)

# Unit test for image_mapping_utils.c
# This unit test depends on trajectory_module.o because image_map calls compute_trajectory_crossings_only
UT_IMGMAP_EXEC = test_image_mapping_utils
UT_IMGMAP_SRC = image_mapping_utils.c
$(UT_IMGMAP_EXEC): $(UT_IMGMAP_SRC) trajectory_module.o $(HEADER)
	@echo "Building unit test $(UT_IMGMAP_EXEC)..."
	$(CC) $(CFLAGS) -DUNIT_TEST_IMAGE_MAPPING_UTILS $(UT_IMGMAP_SRC) trajectory_module.o -o $@ $(LDFLAGS)
	@echo "Running unit test $(UT_IMGMAP_EXEC)..."
	./$(UT_IMGMAP_EXEC)

# Unit test for results_generation.c
# Depends on image_mapping_utils.o, trajectory_module.o, photon_chunk_io.o
UT_RESGEN_EXEC = test_results_generation
UT_RESGEN_SRC = results_generation.c
UT_RESGEN_DEPS_OBJS = image_mapping_utils.o trajectory_module.o photon_chunk_io.o 
$(UT_RESGEN_EXEC): $(UT_RESGEN_SRC) $(UT_RESGEN_DEPS_OBJS) $(HEADER)
	@echo "Building unit test $(UT_RESGEN_EXEC)..."
	$(CC) $(CFLAGS) -DUNIT_TEST_RESULTS_GENERATION $(UT_RESGEN_SRC) $(UT_RESGEN_DEPS_OBJS) -o $@ $(LDFLAGS)
	@echo "Running unit test $(UT_RESGEN_EXEC)..."
	./$(UT_RESGEN_EXEC)

# Unit test for photon_rendering.c (map_photons)
# Depends on ppm_image_io.o, photon_chunk_io.o
UT_PHOTREND_EXEC = test_photon_rendering
UT_PHOTREND_SRC = photon_rendering.c
UT_PHOTREND_DEPS_OBJS = ppm_image_io.o photon_chunk_io.o
$(UT_PHOTREND_EXEC): $(UT_PHOTREND_SRC) $(UT_PHOTREND_DEPS_OBJS) $(HEADER)
	@echo "Building unit test $(UT_PHOTREND_EXEC)..."
	$(CC) $(CFLAGS) -DUNIT_TEST_PHOTON_RENDERING $(UT_PHOTREND_SRC) $(UT_PHOTREND_DEPS_OBJS) -o $@ $(LDFLAGS)
	@echo "Running unit test $(UT_PHOTREND_EXEC)..."
	./$(UT_PHOTREND_EXEC)

# Unit test for photon_chunk_io.c
UT_CHUNKIO_EXEC = test_photon_chunk_io
UT_CHUNKIO_SRC = photon_chunk_io.c
$(UT_CHUNKIO_EXEC): $(UT_CHUNKIO_SRC) $(HEADER)
	@echo "Building unit test $(UT_CHUNKIO_EXEC)..."
	$(CC) $(CFLAGS) -DUNIT_TEST_CHUNK_IO $(UT_CHUNKIO_SRC) -o $@ $(LDFLAGS)
	@echo "Running unit test $(UT_CHUNKIO_EXEC)..."
	./$(UT_CHUNKIO_EXEC)

# Unit test for ppm_image_io.c
UT_PPMIO_EXEC = test_ppm_image_io
UT_PPMIO_SRC = ppm_image_io.c
$(UT_PPMIO_EXEC): $(UT_PPMIO_SRC) $(HEADER)
	@echo "Building unit test $(UT_PPMIO_EXEC)..."
	$(CC) $(CFLAGS) -DUNIT_TEST_PPM_IO $(UT_PPMIO_SRC) -o $@ $(LDFLAGS)
	@echo "Running unit test $(UT_PPMIO_EXEC)..."
	./$(UT_PPMIO_EXEC)

# Target to run all defined unit tests
run_unit_tests: $(UT_TRAJ_EXEC) $(UT_IMGMAP_EXEC) $(UT_RESGEN_EXEC) $(UT_PHOTREND_EXEC) $(UT_CHUNKIO_EXEC) $(UT_PPMIO_EXEC)
	@echo "\nAll specified unit tests have been run."

# Clean target
clean:
	@echo "Cleaning up..."
	rm -f $(OBJS_LIB) $(MAIN_TEST_OBJ) $(MAIN_TEST_EXEC) 
	rm -f $(UT_TRAJ_EXEC) $(UT_IMGMAP_EXEC) $(UT_RESGEN_EXEC) $(UT_PHOTREND_EXEC) $(UT_CHUNKIO_EXEC) $(UT_PPMIO_EXEC)
	rm -f *.o # Catch any other .o files
	rm -f test_*.ppm test_*.bin # Clean up test output files
	rm -f ResCart_*.bin ResRad_*.bin LightRingSeg0thLrg_*.bin # Clean up results_... output files
	rm -f test_lensed_*.ppm # Clean up map_photons output from main_test
	@echo "Cleanup complete."

# Phony targets are not files
.PHONY: all clean run_unit_tests $(UT_TRAJ_EXEC) $(UT_IMGMAP_EXEC) $(UT_RESGEN_EXEC) $(UT_PHOTREND_EXEC) $(UT_CHUNKIO_EXEC) $(UT_PPMIO_EXEC)

## END of Makefile
## Bottom