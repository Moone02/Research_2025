// photon_chunk_io.c
// Top

#define _POSIX_C_SOURCE 200809L
#include "schwarzschild_tracer.h" // For PhotonMapDataPoint, and common includes like stdio, stdlib, stdint
// #include <stdio.h>   // Already in schwarzschild_tracer.h
// #include <stdlib.h>  // Already in schwarzschild_tracer.h
// #include <stdint.h>  // Already in schwarzschild_tracer.h

// --- Chunk I/O (Binary) ---
int save_photon_data_chunk(const char *filename, const PhotonMapDataPoint* data, size_t num_points) {
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        perror("save_photon_data_chunk: fopen failed");
        return -1;
    }
    // Check for num_points overflow before casting to uint32_t
    if (num_points > UINT32_MAX) {
        fprintf(stderr, "save_photon_data_chunk: num_points (%zu) exceeds uint32_t max.\n", num_points);
        fclose(fp);
        return -2; // Indicate different error for overflow
    }
    uint32_t num_points_u32 = (uint32_t)num_points;

    if (fwrite(&num_points_u32, sizeof(uint32_t), 1, fp) != 1) {
        perror("save_photon_data_chunk: fwrite num_points failed");
        fclose(fp);
        return -1;
    }
    if (num_points > 0) { // Only write data if there are points
        if (!data) { // Safety check if num_points > 0 but data is NULL
            fprintf(stderr, "save_photon_data_chunk: num_points > 0 but data is NULL.\n");
            fclose(fp);
            return -3; // Indicate invalid argument
        }
        if (fwrite(data, sizeof(PhotonMapDataPoint), num_points, fp) != num_points) {
            perror("save_photon_data_chunk: fwrite data points failed");
            fclose(fp);
            return -1;
        }
    }
    if (fclose(fp) != 0) {
        perror("save_photon_data_chunk: fclose failed");
        return -1; 
    }
    return 0;
}

PhotonMapDataPoint* load_photon_data_chunk(const char *filename, size_t *num_points_read) {
    if (!num_points_read) { 
        fprintf(stderr, "load_photon_data_chunk: num_points_read pointer is NULL.\n");
        return NULL; 
    }
    *num_points_read = 0; 

    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        perror("load_photon_data_chunk: fopen failed");
        return NULL;
    }
    uint32_t num_points_file;
    if (fread(&num_points_file, sizeof(uint32_t), 1, fp) != 1) {
        // Could be EOF if file is empty or < 4 bytes, not necessarily a full perror.
        if (feof(fp)) {
             fprintf(stderr, "load_photon_data_chunk: Premature EOF reading num_points from %s.\n", filename);
        } else {
            perror("load_photon_data_chunk: fread num_points failed");
        }
        fclose(fp);
        return NULL;
    }
    *num_points_read = (size_t)num_points_file;
    if (*num_points_read == 0) {
        fclose(fp);
        return NULL; 
    }

    PhotonMapDataPoint *data = malloc(*num_points_read * sizeof(PhotonMapDataPoint));
    if (!data) {
        perror("load_photon_data_chunk: malloc failed");
        fclose(fp);
        *num_points_read = 0; // Reset on failure
        return NULL;
    }
    if (fread(data, sizeof(PhotonMapDataPoint), *num_points_read, fp) != *num_points_read) {
        perror("load_photon_data_chunk: fread data points failed or short read");
        free(data);
        fclose(fp);
        *num_points_read = 0; // Reset on failure
        return NULL;
    }
    if (fclose(fp) != 0) {
        perror("load_photon_data_chunk: fclose failed");
    }
    return data;
}


#ifdef UNIT_TEST_CHUNK_IO
#include <math.h> // For fabs in comparison helper

// Helper to compare PhotonMapDataPoint for testing
static int test_compare_photons(const PhotonMapDataPoint* p1, const PhotonMapDataPoint* p2, double tol) {
    return (fabs(p1->y_window_cart_x - p2->y_window_cart_x) < tol &&
            fabs(p1->y_window_cart_y - p2->y_window_cart_y) < tol &&
            fabs(p1->y_image_cart_x - p2->y_image_cart_x) < tol &&
            fabs(p1->y_image_cart_y - p2->y_image_cart_y) < tol);
}

int main() {
    printf("--- Unit Test for photon_chunk_io ---\n");
    const char* test_file = "test_chunk_io_temp.bin";
    PhotonMapDataPoint *loaded_data;
    size_t num_loaded;

    // Test 1: Save and load 3 points
    PhotonMapDataPoint original_data[3] = {
        {1.0, 2.0, 3.0, 4.0},
        {-1.0, -2.0, -3.0, -4.0},
        {0.5, -0.5, 10.5, -10.5}
    };
    printf("Test 1: Saving 3 points...\n");
    if (save_photon_data_chunk(test_file, original_data, 3) != 0) {
        printf("Test 1 FAILED: save_photon_data_chunk\n"); return 1;
    }
    loaded_data = load_photon_data_chunk(test_file, &num_loaded);
    if (!loaded_data || num_loaded != 3) {
        printf("Test 1 FAILED: load_photon_data_chunk (count %zu, ptr %p)\n", num_loaded, (void*)loaded_data);
        if(loaded_data) free(loaded_data); remove(test_file); return 1;
    }
    int match = 1;
    for(size_t i=0; i<3; ++i) if(!test_compare_photons(&original_data[i], &loaded_data[i], 1e-9)) match=0;
    if(match) printf("Test 1 PASSED\n"); else printf("Test 1 FAILED: data mismatch\n");
    free(loaded_data);
    remove(test_file);

    // Test 2: Save and load 0 points
    printf("Test 2: Saving 0 points...\n");
    if (save_photon_data_chunk(test_file, NULL, 0) != 0) { // Pass NULL for data when 0 points
        printf("Test 2 FAILED: save_photon_data_chunk (0 points)\n"); return 1;
    }
    loaded_data = load_photon_data_chunk(test_file, &num_loaded);
    if (loaded_data == NULL && num_loaded == 0) {
        printf("Test 2 PASSED (0 points handling)\n");
    } else {
        printf("Test 2 FAILED: load_photon_data_chunk (0 points) (count %zu, ptr %p)\n", num_loaded, (void*)loaded_data);
        if(loaded_data) free(loaded_data);
    }
    remove(test_file);
    
    // Test 3: Attempt to load non-existent file
    printf("Test 3: Loading non-existent file...\n");
    loaded_data = load_photon_data_chunk("non_existent_chunk.bin", &num_loaded);
    if (loaded_data == NULL && num_loaded == 0) {
         printf("Test 3 PASSED (non-existent file handling)\n");
    } else {
         printf("Test 3 FAILED: load non-existent file did not return NULL/0 (count %zu, ptr %p)\n", num_loaded, (void*)loaded_data);
        if(loaded_data) free(loaded_data);
    }


    printf("--- Unit Test for photon_chunk_io Finished ---\n");
    return 0;
}
#endif // UNIT_TEST_CHUNK_IO

// photon_chunk_io.c
// Bottom