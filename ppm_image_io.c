// ppm_image_io.c
// Top

#define _POSIX_C_SOURCE 200809L // For strdup, though not directly used here, good practice if other utils need it
#include "schwarzschild_tracer.h" // For PPMImage, RGBColor, and other common includes like stdio, stdlib, ctype, string
// #include <stdio.h>  // Already in schwarzschild_tracer.h
// #include <stdlib.h> // Already in schwarzschild_tracer.h
// #include <string.h> // Already in schwarzschild_tracer.h
// #include <ctype.h>  // Already in schwarzschild_tracer.h

// --- PPM Image I/O (Basic P6 support) ---
static void skip_ppm_comments_and_whitespace(FILE *fp) {
    int ch;
    // Skip leading whitespace
    while ((ch = fgetc(fp)) != EOF && isspace(ch)) {
        // Keep reading
    }
    // If it's a comment, skip to end of line
    if (ch == '#') {
        while ((ch = fgetc(fp)) != EOF && ch != '\n' && ch != '\r') {
            // Keep reading comment
        }
        // After skipping comment, there might be more whitespace or another comment
        skip_ppm_comments_and_whitespace(fp);
    } else if (ch != EOF) {
        // If it wasn't whitespace or '#', put it back
        ungetc(ch, fp);
    }
}

PPMImage* load_ppm_image(const char *filename) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        perror("load_ppm_image: fopen failed");
        return NULL;
    }

    char magic[3];
    if (fgets(magic, sizeof(magic), fp) == NULL || strncmp(magic, "P6", 2) != 0) {
        fprintf(stderr, "load_ppm_image: Not a P6 PPM file or read error (magic: %s).\n", magic);
        fclose(fp);
        return NULL;
    }

    PPMImage *image = calloc(1, sizeof(PPMImage));
    if (!image) {
        perror("load_ppm_image: calloc PPMImage failed");
        fclose(fp);
        return NULL;
    }

    skip_ppm_comments_and_whitespace(fp);
    if (fscanf(fp, "%d", &image->width) != 1) {
        fprintf(stderr, "load_ppm_image: Failed to read width.\n");
        goto read_error_ppm_load;
    }
    skip_ppm_comments_and_whitespace(fp);
    if (fscanf(fp, "%d", &image->height) != 1) {
        fprintf(stderr, "load_ppm_image: Failed to read height.\n");
        goto read_error_ppm_load;
    }
    skip_ppm_comments_and_whitespace(fp);
    int max_val;
    if (fscanf(fp, "%d", &max_val) != 1 || max_val != 255) {
        fprintf(stderr, "load_ppm_image: PPM max color value not 255 or read error (max_val: %d).\n", max_val);
        goto read_error_ppm_load;
    }

    int ch_after_maxval = fgetc(fp);
    if (!isspace(ch_after_maxval)) {
        fprintf(stderr, "load_ppm_image: Expected whitespace after max color value, got '%c' (ASCII: %d).\n", ch_after_maxval, ch_after_maxval);
        goto read_error_ppm_load;
    }

    if (image->width <= 0 || image->height <= 0) {
        fprintf(stderr, "load_ppm_image: Invalid image dimensions %dx%d.\n", image->width, image->height);
        goto read_error_ppm_load;
    }

    image->channels = 3;
    size_t data_size = (size_t)image->width * image->height * image->channels;
    if (data_size == 0) { // Handle 0xN or Nx0 images
        image->data = NULL;
    } else {
        image->data = malloc(data_size);
        if (!image->data) {
            perror("load_ppm_image: malloc for pixel data failed");
            goto read_error_ppm_load;
        }
        if (fread(image->data, sizeof(unsigned char), data_size, fp) != data_size) {
            perror("load_ppm_image: fread for pixel data failed or short read");
            free(image->data); image->data = NULL; 
            goto read_error_ppm_load;
        }
    }

    if (fclose(fp) != 0) {
        perror("load_ppm_image: fclose failed");
    }
    return image;

read_error_ppm_load:
    if (fp) { fclose(fp); }
    free_ppm_image(image); 
    return NULL;
}

int save_ppm_image(const char *filename, const PPMImage *image) {
    if (!image || image->width <= 0 || image->height <= 0 || image->channels != 3) {
        // Allow image->data to be NULL if width or height is 0
        if (!(image && (image->width == 0 || image->height == 0) && image->data == NULL)) {
             fprintf(stderr, "save_ppm_image: Invalid image data for saving.\n");
             return -1;
        }
    }
    // If dimensions are zero, but data is not NULL (or vice-versa, excluding valid 0-dim case)
    if ((image->width == 0 || image->height == 0) && image->data != NULL) {
        fprintf(stderr, "save_ppm_image: Zero dimensions but non-NULL data.\n");
        return -1;
    }
    if ((image->width > 0 && image->height > 0) && image->data == NULL) {
        fprintf(stderr, "save_ppm_image: Non-zero dimensions but NULL data.\n");
        return -1;
    }


    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        perror("save_ppm_image: fopen failed");
        return -1;
    }

    if (fprintf(fp, "P6\n%d %d\n255\n", image->width, image->height) < 0) {
        perror("save_ppm_image: fprintf header failed");
        fclose(fp);
        return -1;
    }

    if (image->width > 0 && image->height > 0 && image->data) { // Only write data if there is any
        size_t data_size = (size_t)image->width * image->height * image->channels;
        if (fwrite(image->data, sizeof(unsigned char), data_size, fp) != data_size) {
            perror("save_ppm_image: fwrite pixel data failed");
            fclose(fp);
            return -1;
        }
    }

    if (fclose(fp) != 0) {
        perror("save_ppm_image: fclose failed");
        return -1; 
    }
    return 0;
}

void free_ppm_image(PPMImage *image) {
    if (!image) { return; }
    if(image->data) { free(image->data); image->data = NULL; } // Good practice to NULL after free
    free(image);
}


#ifdef UNIT_TEST_PPM_IO
// Helper to create a dummy PPM for testing
static PPMImage* test_create_dummy_ppm(int width, int height, unsigned char val) {
    PPMImage* img = malloc(sizeof(PPMImage));
    img->width = width; img->height = height; img->channels = 3;
    if (width == 0 || height == 0) {
        img->data = NULL;
        return img;
    }
    img->data = malloc((size_t)width * height * 3);
    for(int i=0; i < width*height*3; ++i) img->data[i] = val + (i % 50); // Some variation
    return img;
}

// Compare PPMImages (simplified for test)
static int test_compare_ppm(const PPMImage* img1, const PPMImage* img2) {
    if (!img1 && !img2) return 1;
    if (!img1 || !img2) return 0;
    if (img1->width != img2->width || img1->height != img2->height || img1->channels != img2->channels) return 0;
    if ((img1->data == NULL) != (img2->data == NULL)) return 0;
    if (img1->data) { // If data is not NULL, compare it
        return memcmp(img1->data, img2->data, (size_t)img1->width * img1->height * img1->channels) == 0;
    }
    return 1; // Both NULL data, dimensions match
}


int main() {
    printf("--- Unit Test for ppm_image_io ---\n");
    const char* test_file = "test_ppm_io_temp.ppm";
    PPMImage *img_orig, *img_loaded;

    // Test 1: Save and load a small image
    img_orig = test_create_dummy_ppm(5, 3, 100);
    if (!img_orig) { printf("Test 1 FAILED: dummy creation failed\n"); return 1;}
    printf("Test 1: Saving 5x3 image...\n");
    if (save_ppm_image(test_file, img_orig) != 0) {
        printf("Test 1 FAILED: save_ppm_image\n");
        free_ppm_image(img_orig); return 1;
    }
    img_loaded = load_ppm_image(test_file);
    if (!img_loaded) {
        printf("Test 1 FAILED: load_ppm_image\n");
        free_ppm_image(img_orig); remove(test_file); return 1;
    }
    if (test_compare_ppm(img_orig, img_loaded)) {
        printf("Test 1 PASSED\n");
    } else {
        printf("Test 1 FAILED: images differ\n");
    }
    free_ppm_image(img_orig);
    free_ppm_image(img_loaded);
    remove(test_file);

    // Test 2: Zero dimension image
    img_orig = test_create_dummy_ppm(0, 10, 0); // 0 width
    printf("Test 2: Saving 0x10 image...\n");
    if (save_ppm_image(test_file, img_orig) != 0) {
        printf("Test 2 FAILED: save_ppm_image (0-width)\n");
        free_ppm_image(img_orig); return 1;
    }
    img_loaded = load_ppm_image(test_file);
     if (!img_loaded) {
        printf("Test 2 FAILED: load_ppm_image (0-width) returned NULL\n");
        free_ppm_image(img_orig); remove(test_file); return 1;
    }
    if (img_loaded->width == 0 && img_loaded->height == 10 && img_loaded->data == NULL) {
        printf("Test 2 PASSED (0-width image handling)\n");
    } else {
        printf("Test 2 FAILED: 0-width image loaded incorrectly (w:%d,h:%d,data:%p)\n", img_loaded->width, img_loaded->height, (void*)img_loaded->data);
    }
    free_ppm_image(img_orig);
    free_ppm_image(img_loaded);
    remove(test_file);

    printf("--- Unit Test for ppm_image_io Finished ---\n");
    return 0;
}
#endif // UNIT_TEST_PPM_IO

// ppm_image_io.c
// Bottom