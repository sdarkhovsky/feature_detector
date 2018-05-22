// feature_detector.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "lpngwrapper.hpp"
#include <list>
#include <vector>

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

void calculate_histogram(MatrixXd& rgb_channel, int x0, int y0, int wx, int wy, VectorXd& hist, int num_bins)
{
    int bin;
    int x, y;

    hist = VectorXd::Zero(num_bins);
    int width = rgb_channel.cols();
    int height = rgb_channel.rows();
    double bin_size = 1.0 / num_bins;

    if (y0 < wy || y0 + wy > height || x0 < wx || x0 + wx > width)
        return;

    for (y = y0 - wy; y < y0+wy; y++) {
        for (x = x0-wx; x < x0+wx; x++) {
            bin = rgb_channel(y, x)/ bin_size;
            if (bin >= num_bins)
                bin = num_bins - 1;
            hist(bin)++;
        }
    }
    hist /= ((2 * wx + 1)*(2 * wy + 1));
}

void compare_histograms(VectorXd histogram_1[], std::vector<std::vector<VectorXd>> histograms_2[], int num_channels, std::string vis_image_path)
{
    int i;
    int x, y;
    double min_val, max_val;
    VectorXd hist, hist1, hist2;
    MatrixXd vis_rgb_channels[3];
    int height = histograms_2[0].size();
    int width = histograms_2[0][0].size();
    double val;

    for (i = 0; i < num_channels; i++)
    {
        vis_rgb_channels[i] = MatrixXd::Zero(height, width);
    }

    for (i = 0; i < num_channels; i++)
    {
        for (y = 0; y < height; y++)
        {
            for (x = 0; x < width; x++)
            {
                hist = histogram_1[i] - histograms_2[i][y][x];

                //1111111111111111111111111111111111
                max_val = histogram_1[i].maxCoeff();
                min_val = histogram_1[i].minCoeff();

                max_val = histograms_2[i][y][x].maxCoeff();
                min_val = histograms_2[i][y][x].minCoeff();

                max_val = hist.maxCoeff();
                min_val = hist.minCoeff();

                if (x == 75 && y == 37)
                {
                    val = 0;
                }
                //11111111111111111111111111111111111

                val = hist.norm();
                val /= hist.size();
                val *= 255.0;
                vis_rgb_channels[i](y, x) = val;
            }
        }
    }

    errno_t err = write_png_file(vis_image_path.c_str(), vis_rgb_channels, 3);
}

int main(int argc, char **argv)
{
    MatrixXd rgb_channels_1[3], rgb_channels_2[3], rgb_channels_3[3];
    VectorXd histogram_1[3];
    std::vector<std::vector<VectorXd>> histograms_2[3];
    int num_channels = 3;

    std::string image_path[3];
    std::string arg;
    int x0, y0, wx, wy;
    int num_bins;

/*
    Command line arguments:
-num_bins 10 -wx 3 -wy 3 -x0 75 -y0 37 -i1 C:\Users\Sam\Documents\EyeSignalsProjects\images\image1_025.png -i2 C:\Users\Sam\Documents\EyeSignalsProjects\images\image2_025.png -o1 C:\Users\Sam\Documents\EyeSignalsProjects\images\image_out.png
*/
    int i = 1;
    while(i < argc)
    {
        arg = argv[i];
        if (arg == "-i1")
        {
            i++;
            image_path[0] = argv[i];
        }
        if (arg == "-i2")
        {
            i++;
            image_path[1] = argv[i];
        }
        if (arg == "-o1")
        {
            i++;
            image_path[2] = argv[i];
        }

        if (arg == "-x0")
        {
            i++;
            x0 = atoi(argv[i]);
        }

        if (arg == "-y0")
        {
            i++;
            y0 = atoi(argv[i]);
        }

        if (arg == "-wx")
        {
            i++;
            wx = atoi(argv[i]);
        }

        if (arg == "-wy")
        {
            i++;
            wy = atoi(argv[i]);
        }

        if (arg == "-num_bins")
        {
            i++;
            num_bins = atoi(argv[i]);
        }
        i++;
    }

    read_png_file(image_path[0].c_str(), rgb_channels_1, num_channels);
    read_png_file(image_path[1].c_str(), rgb_channels_2, num_channels);
    prepare_channels(rgb_channels_1, num_channels);
    prepare_channels(rgb_channels_2, num_channels);

    if (x0 < 0 || x0 >= rgb_channels_1[0].cols())
    {
        fprintf(stderr, "x0 out of range: %d\n", x0);
        exit(1);
    }

    if (y0 < 0 || y0 >= rgb_channels_1[0].rows())
    {
        fprintf(stderr, "y0 out of range: %d\n", y0);
        exit(1);
    }

    for (i = 0; i < 3; i++)
    {
        calculate_histogram(rgb_channels_1[i], x0, y0, wx, wy, histogram_1[i], num_bins);
    }

    int width = rgb_channels_2[0].cols();
    int height = rgb_channels_2[0].rows();
    int x, y;

    for (i = 0; i < 3; i++)
    {
        histograms_2[i].resize(height);
        for (y = 0; y < height; y++) {
            histograms_2[i][y].resize(width);
            for (x = 0; x < width; x++) {
                calculate_histogram(rgb_channels_2[i], x, y, wx, wy, histograms_2[i][y][x], num_bins);
            }
        }
    }

    compare_histograms(histogram_1, histograms_2, num_channels, image_path[2]);
//    MatrixXd image_diff;
//    compare_channels(rgb_channels_1, rgb_channels_2, num_channels, image_diff, image_path[2]);

    return 0;
}

