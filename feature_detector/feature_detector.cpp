// feature_detector.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "lpngwrapper.hpp"
#include <list>


void usage(void)
{
    fprintf(stderr, "usage: feature_detector <input_image_1_file_path> <input_image_2_file_path> <output_image_3_file_path>\n");
    exit(1);
}

void prepare_channels(MatrixXd rgb_channels[], int num_channels)
{
    int x, y;
    double r,g,b,sum;

    int width = rgb_channels[0].cols();
    int height = rgb_channels[0].rows();

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            r = rgb_channels[0](y, x);
            g = rgb_channels[1](y, x);
            b = rgb_channels[2](y, x);
            sum = r + g + b;
            rgb_channels[0](y, x) = r / sum;
            rgb_channels[1](y, x) = g / sum;
            rgb_channels[2](y, x) = b / sum;
        }
    }
}

void compare_channels(MatrixXd rgb_channels1[], MatrixXd rgb_channels2[], int num_channels, MatrixXd& image_diff, std::string vis_image_path)
{
    int x1, y1, x2, y2;
    int sh_x2, sh_y2;  // shift of the second image
    int sh_x, sh_y;
    int i;
    int diff_cnt;
    double ch_diff, sum;
    int width = rgb_channels1[0].cols();
    int height = rgb_channels1[0].rows();

    int sh_step_x = 4, sh_step_y = 4;
    int sh_width = width/sh_step_x;
    int sh_height = height/sh_step_y;

    image_diff = MatrixXd::Zero(sh_height, sh_width);

    std::list<std::pair <int, int>> min_diff_list;

    /* compare the shifted image 2 with the image 1*/
    for (sh_y = 0; sh_y < sh_height; sh_y++) {
        sh_y2 = (sh_y - sh_height / 2)*sh_step_y;
        for (sh_x = 0; sh_x < sh_width; sh_x++) {
            sh_x2 = (sh_x - sh_width / 2)*sh_step_x;
            /* element (x1,y1) of the image 1 is compared against (x1-sh_x2,y1-sh_y2) of the image 2 */
            sum = 0;
            diff_cnt = 0;
            for (y1 = 0; y1 < height; y1++) {
                y2 = y1 - sh_y2;
                for (x1 = 0; x1 < width; x1++) {
                    x2 = x1 - sh_x2;
                    if (x2 < 0 || x2 >= width || y2 < 0 || y2 >= height)
                        continue;
                    for (i = 0; i < 3; i++) {
                        ch_diff = rgb_channels1[i](y1, x1) - rgb_channels2[i](y2, x2);
                        sum += ch_diff * ch_diff;
                        diff_cnt++;
                    }
                }
            }

            if (diff_cnt > 0)
                image_diff(sh_y, sh_x) = sqrt(sum) / diff_cnt;
            else
                image_diff(sh_y, sh_x) = FLT_MAX;
        }
    }

    // calculate local min in image_diff
    for (sh_y = 1; sh_y < sh_height-1; sh_y++) {
//        printf("\n");
        for (sh_x = 1; sh_x < sh_width-1; sh_x++) {
//            printf("%f ", image_diff(sh_y, sh_x));

            if (image_diff(sh_y, sh_x) < image_diff(sh_y-1, sh_x) &&
                image_diff(sh_y, sh_x) < image_diff(sh_y, sh_x-1) &&
                image_diff(sh_y, sh_x) < image_diff(sh_y+1, sh_x) &&
                image_diff(sh_y, sh_x) < image_diff(sh_y, sh_x + 1))
            {
                printf("sh_x=%d sh_y=%d\n", sh_x, sh_y);
                min_diff_list.push_back(std::make_pair(sh_x, sh_y));
            }
        }
    }

    assert(num_channels == 3);
    MatrixXd vis_rgb_channels[3];
    for (i = 0; i < num_channels; i++)
    {
        vis_rgb_channels[i] = MatrixXd::Zero(height, width);
    }

    std::size_t found = vis_image_path.find_last_of(".");
    assert(found != std::string::npos);
    std::string vis_image_path_base = vis_image_path.substr(0, found);
    for (std::list<std::pair <int, int>>::iterator it = min_diff_list.begin(); it != min_diff_list.end(); ++it)
    {
        sh_x = it->first;
        sh_y = it->second;
        sh_y2 = (sh_y - sh_height / 2)*sh_step_y;
        sh_x2 = (sh_x - sh_width / 2)*sh_step_x;
        for (y1 = 0; y1 < height; y1++) {
            y2 = y1 - sh_y2;
            for (x1 = 0; x1 < width; x1++) {
                x2 = x1 - sh_x2;
                if (x2 < 0 || x2 >= width || y2 < 0 || y2 >= height)
                {
                    for (i = 0; i < 3; i++)
                        vis_rgb_channels[i](y1, x1) = 255.0;
                }
                else
                {
                    sum = 0;
                    for (i = 0; i < 3; i++) {
                        ch_diff = rgb_channels1[i](y1, x1) - rgb_channels2[i](y2, x2);
                        sum += ch_diff * ch_diff;
                    }
                    for (i = 0; i < 3; i++)
                        vis_rgb_channels[i](y1, x1) = sum*255.0;
                }
            }
        }

        std::string image_path = vis_image_path_base + "_" + std::to_string(sh_x) + "_" + std::to_string(sh_y) + ".png";
        errno_t err = write_png_file(image_path.c_str(), vis_rgb_channels, 3);
    }
}

int main(int argc, char **argv)
{
    MatrixXd rgb_channels_1[3], rgb_channels_2[3], rgb_channels_3[3];
    int num_channels = 3;

    if (argc != 4)
    {
        usage();
    }

    std::string image1_path = argv[1];
    std::string image2_path = argv[2];
    std::string image3_path = argv[3];

    read_png_file(image1_path.c_str(), rgb_channels_1, num_channels);
    read_png_file(image2_path.c_str(), rgb_channels_2, num_channels);
    prepare_channels(rgb_channels_1, num_channels);
    prepare_channels(rgb_channels_2, num_channels);

    MatrixXd image_diff;
    compare_channels(rgb_channels_1, rgb_channels_2, num_channels, image_diff, image3_path);

    return 0;
}

