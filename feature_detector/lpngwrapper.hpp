#ifndef LPNGWRAPPER_HPP
#define LPNGWRAPPER_HPP



#include <stdio.h>
#include <assert.h>
#include <png.h>
#include <string.h>
#include <stdlib.h>

#include <Eigen/Dense>

using namespace Eigen;

typedef struct {
    int width;
    int height;
    png_size_t row_size;
    png_size_t size;
    int color_type;
    png_byte* data;
} RawImageData;
 
typedef struct {
    png_byte* data;
    png_size_t size;
    png_size_t offset;
} ReadDataHandle;
 
typedef struct {
    png_uint_32 width;
    png_uint_32 height;
    int color_type;
} PngInfo;

class CPNGWrapper {
public:
    CPNGWrapper();
    ~CPNGWrapper();

    void get_png_channels(MatrixXd rgb_channels[], int num_channels);

    /* write functions */
    bool png_write_prepare(FILE *fp, int width, int height, png_structp* png_ptrptr, png_infop* info_ptrptr);
    bool png_write_complete(png_structp png_ptr, png_bytep * row_pointers);

    /* Returns the decoded image data, or aborts if there's an error during decoding. */
    void get_raw_image_data_from_png();
    PngInfo read_and_update_info(const png_structp png_ptr, const png_infop info_ptr);
    void read_entire_png_image(const png_structp png_ptr, const png_infop info_ptr, const png_uint_32 height);
    RawImageData m_raw_image_data;
    void* m_pMemoryImage;
    size_t m_iMemoryImageLen;

    bool m_bReleaseFileMemory;
};

errno_t read_png_file(const char* file_path, MatrixXd rgb_channels[], int num_channels);
errno_t write_png_file(const char* file_name, MatrixXd rgb_channels[], int num_channels);

#endif // LPNGWRAPPER_HPP
