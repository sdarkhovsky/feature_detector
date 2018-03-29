// feature_detector.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "lpngwrapper.hpp"

void usage(void)
{
    fprintf(stderr, "usage: feature_detector <input_image_file_path> <output_image_file_path>\n");
    exit(1);
}

int main(int argc, char **argv)
{
    MatrixXd rgb_channels[3];

    if (argc != 3)
    {
        usage();
    }

    read_png_file(argv[1], rgb_channels, 3);
    rgb_channels[0].setZero();
    rgb_channels[2].setZero();
    errno_t err = write_png_file(argv[2], rgb_channels, 3);

    return 0;
}

