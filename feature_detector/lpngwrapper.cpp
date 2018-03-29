#include "stdafx.h"

#include "lpngwrapper.hpp"
#include <vector>
#include <memory>

// should the png file's color type be RGBA 8 bits per channel?
static void read_png_data_callback(
    png_structp png_ptr, png_byte* raw_data, png_size_t read_length) {
    ReadDataHandle* handle = (ReadDataHandle*)png_get_io_ptr(png_ptr);
    const png_byte* png_src = handle->data + handle->offset;
 
    memcpy(raw_data, png_src, read_length);
    handle->offset += read_length;
}

CPNGWrapper::CPNGWrapper()
{
    m_bReleaseFileMemory = false;

    m_raw_image_data.data = nullptr;
}

CPNGWrapper::~CPNGWrapper()
{
    if (m_bReleaseFileMemory)
        free(m_pMemoryImage);

    //todo: the data can be released earlier: after they were converted to YUVA
    if (m_raw_image_data.data)
        free((void*)m_raw_image_data.data);
}

PngInfo CPNGWrapper::read_and_update_info(const png_structp png_ptr, const png_infop info_ptr)
{
    int bit_depth, color_type;
    PngInfo pngInfo;

    png_read_info(png_ptr, info_ptr);
    png_get_IHDR(
        png_ptr, info_ptr, &pngInfo.width, &pngInfo.height, &bit_depth, &color_type, NULL, NULL, NULL);

    // Convert transparency to full alpha
    if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS))
        png_set_tRNS_to_alpha(png_ptr);

    // Convert grayscale, if needed.
    if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
        png_set_expand_gray_1_2_4_to_8(png_ptr);

    // Convert paletted images, if needed.
    if (color_type == PNG_COLOR_TYPE_PALETTE)
        png_set_palette_to_rgb(png_ptr);

    // Add alpha channel, if there is none.
    // Rationale: GL_RGBA is faster than GL_RGB on many GPUs)
    if (color_type == PNG_COLOR_TYPE_PALETTE || color_type == PNG_COLOR_TYPE_RGB)
        png_set_add_alpha(png_ptr, 0xFF, PNG_FILLER_AFTER);

    // Ensure 8-bit packing
    if (bit_depth < 8)
        png_set_packing(png_ptr);
    /* todo: temporarily commented out because png_set_scale_16 is not defined
        else if (bit_depth == 16)
            png_set_scale_16(png_ptr);
    */

    png_read_update_info(png_ptr, info_ptr);

    // Read the new color type after updates have been made.
    pngInfo.color_type = png_get_color_type(png_ptr, info_ptr);

    return pngInfo;
}

void CPNGWrapper::read_entire_png_image(
    const png_structp png_ptr,
    const png_infop info_ptr,
    const png_uint_32 height)
{
    const png_size_t row_size = png_get_rowbytes(png_ptr, info_ptr);
    const png_size_t data_length = row_size * height;
    assert(row_size > 0);

    png_byte* raw_image = (png_byte*)malloc(data_length);
    if (raw_image != NULL)
    {
        std::vector <png_byte*> row_ptrs;
        row_ptrs.resize(height);

        png_uint_32 i;
        for (i = 0; i < height; i++) {
            row_ptrs[i] = raw_image + i * row_size;
        }

        png_read_image(png_ptr, &row_ptrs[0]);
    }

    m_raw_image_data.row_size = row_size;
    m_raw_image_data.size = data_length;
    m_raw_image_data.data = raw_image;
}


void CPNGWrapper::get_raw_image_data_from_png()
{
    assert(m_pMemoryImage != NULL && m_iMemoryImageLen > 8);
    assert(png_check_sig((png_bytep)m_pMemoryImage, 8));

    png_structp png_ptr = png_create_read_struct(
        PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    assert(png_ptr != NULL);
    png_infop info_ptr = png_create_info_struct(png_ptr);
    assert(info_ptr != NULL);
    ReadDataHandle png_data_handle;
    png_data_handle.data = (png_byte*)m_pMemoryImage;
    png_data_handle.size = m_iMemoryImageLen;
    png_data_handle.offset = 0;
    png_set_read_fn(png_ptr, &png_data_handle, read_png_data_callback);

    if (setjmp(png_jmpbuf(png_ptr))) {
        assert(false);
    }

    const PngInfo png_info = read_and_update_info(png_ptr, info_ptr);

    read_entire_png_image(png_ptr, info_ptr, png_info.height);

    png_read_end(png_ptr, info_ptr);
    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);

    m_raw_image_data.width = png_info.width;
    m_raw_image_data.height = png_info.height;
    m_raw_image_data.color_type = png_info.color_type;
}

void CPNGWrapper::get_png_channels(MatrixXd rgb_channels[], int num_channels)
{
    assert(m_raw_image_data.color_type == PNG_COLOR_TYPE_RGB_ALPHA);
    assert(num_channels == 3);

    double R, G, B, A;
    int i, j, ch;

    for (ch = 0; ch < num_channels; ch++)
    {
        rgb_channels[ch] = MatrixXd::Zero(m_raw_image_data.height, m_raw_image_data.width);
    }

    for (i = 0; i < m_raw_image_data.height; i++) {
        png_byte* row_ptr = m_raw_image_data.data + i * m_raw_image_data.row_size;
        for (j = 0; j < m_raw_image_data.width; j++) {
            png_byte* pix_ptr = row_ptr + j * 4;
            R = *pix_ptr++;
            G = *pix_ptr++;
            B = *pix_ptr++;
            A = *pix_ptr;
            rgb_channels[0](i, j) = R;
            rgb_channels[1](i, j) = G;
            rgb_channels[2](i, j) = B;
        }
    }
}

bool CPNGWrapper::png_write_prepare(FILE *fp, int width, int height, png_structp* png_ptrptr, png_infop* info_ptrptr) {

    png_byte color_type = PNG_COLOR_TYPE_RGBA;
    png_byte bit_depth = 8;

    png_structp png_ptr;
    png_infop info_ptr;

    /* initialize stuff */
    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

    if (!png_ptr)
        return false;

    info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr)
        return false;

    if (setjmp(png_jmpbuf(png_ptr)))
        return false;

    png_init_io(png_ptr, fp);

    /* write header */
    if (setjmp(png_jmpbuf(png_ptr)))
        return false;

    png_set_IHDR(png_ptr, info_ptr, width, height,
        bit_depth, color_type, PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    png_write_info(png_ptr, info_ptr);

    /* write bytes */
    if (setjmp(png_jmpbuf(png_ptr)))
        return false;

    *png_ptrptr = png_ptr;
    *info_ptrptr = info_ptr;

    return true;
}


bool CPNGWrapper::png_write_complete(png_structp png_ptr, png_bytep * row_pointers) {

    png_write_image(png_ptr, row_pointers);

    /* end write */
    if (setjmp(png_jmpbuf(png_ptr)))
        return false;

    png_write_end(png_ptr, NULL);

    return true;
}


errno_t write_png_file(const char* file_name, MatrixXd rgb_channels[], int num_channels)
{
    int x, y;
    png_structp png_ptr;
    png_infop info_ptr;
    bool result = true;
    png_byte* ptr;
    FILE *fp = NULL;
    png_bytep * row_pointers = NULL;;

    assert(num_channels == 3);

    #define GENERAL_ERROR 1

    CPNGWrapper png_wrapper;

    int width = rgb_channels[0].cols();
    int height = rgb_channels[0].rows();

    /* create file */
    errno_t err = fopen_s(&fp, file_name, "wb");
    if (err)
        return err;

    if (!png_wrapper.png_write_prepare(fp, width, height, &png_ptr, &info_ptr))
        return GENERAL_ERROR;

    row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * height);
    if (!row_pointers)
        goto Cleanup;
    for (y = 0; y < height; y++) {
        row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png_ptr, info_ptr));
        if (!row_pointers[y])
            goto Cleanup;
    }

    for (y = 0; y < height; y++) {
        png_byte* row = row_pointers[y];
        for (x = 0; x < width; x++) {
            ptr = &(row[x * 4]);

            ptr[0] = rgb_channels[0](y,x);
            ptr[1] = rgb_channels[1](y, x);
            ptr[2] = rgb_channels[2](y, x);
            ptr[3] = 255;
        }
    }

    png_wrapper.png_write_complete(png_ptr, row_pointers);

Cleanup:
    /* cleanup heap allocation */
    if (row_pointers) {
        for (y = 0; y < height; y++) {
            if (row_pointers[y])
                free(row_pointers[y]);
        }
        free(row_pointers);
    }

    fclose(fp);

    return err;
}

errno_t read_png_file(const char* file_path, MatrixXd rgb_channels[], int num_channels)
{
    FILE *f = NULL;
    assert(num_channels == 3);
    errno_t err = fopen_s(&f, file_path, "rb");
    if (err == 0)
    {
        CPNGWrapper png_wrapper;
        void * pBuf;
        int bufLen;

        fseek(f, 0, SEEK_END);
        bufLen = ftell(f);
        fseek(f, 0, SEEK_SET);

        pBuf = malloc(bufLen);
        if (pBuf)
        {
            size_t count = 1;
            size_t nRead = fread(pBuf, bufLen, count, f);
            if (nRead == count)
            {
                png_wrapper.m_pMemoryImage = pBuf;
                png_wrapper.m_iMemoryImageLen = bufLen;
                png_wrapper.m_bReleaseFileMemory = true;
            }
            else
            {
                free(pBuf);
            }
        }
        fclose(f);

        png_wrapper.get_raw_image_data_from_png();
        png_wrapper.get_png_channels(rgb_channels, num_channels);

    }

    return err;
}
