// feature_detector.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "lpngwrapper.hpp"
#include <list>
#include <vector>
#include <iostream>

void usage(void)
{
    fprintf(stderr, "usage: feature_detector <input_image_1_file_path> <input_image_2_file_path> <output_image_3_file_path>\n");
    exit(1);
}

void preprocess_channels(MatrixXd rgb_channels[], int num_channels, MatrixXd& intensity)
{
    int x, y;
    double r,g,b,sum;

    int width = rgb_channels[0].cols();
    int height = rgb_channels[0].rows();
    intensity = MatrixXd::Zero(height, width);

    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            r = rgb_channels[0](y, x);
            g = rgb_channels[1](y, x);
            b = rgb_channels[2](y, x);
            sum = r + g + b;
            assert(r >= 0 && g >= 0 && b >= 0);
            if (sum > 0)
            {
                rgb_channels[0](y, x) = r / sum;
                rgb_channels[1](y, x) = g / sum;
                rgb_channels[2](y, x) = b / sum;
            }
            intensity(y, x) = sum/3.0;
        }
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

    if (y0 < wy || y0 + wy >= height || x0 < wx || x0 + wx >= width)
        return;

    for (y = y0 - wy; y <= y0+wy; y++) {
        for (x = x0-wx; x <= x0+wx; x++) {
            bin = rgb_channel(y, x)/ bin_size;

            if (bin >= num_bins)
                bin = num_bins - 1;
            hist(bin)++;
        }
    }

    hist.normalize();

}

void compare_histograms(VectorXd histogram_1[], std::vector<std::vector<VectorXd>> histograms_2[], int num_channels, std::string vis_image_path, MatrixXd& identity_score)
{
    int i;
    int x, y;
    double dist, max_dist;
    VectorXd hist1, hist2;
    MatrixXd vis_rgb_channels[3];
    int height = histograms_2[0].size();
    int width = histograms_2[0][0].size();

    for (i = 0; i < num_channels; i++)
    {
        vis_rgb_channels[i] = MatrixXd::Zero(height, width);
    }

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            max_dist = 0.0;
            for (i = 0; i < num_channels; i++)
            {
                dist = (histogram_1[i] - histograms_2[i][y][x]).cwiseAbs().maxCoeff();
                if (max_dist < dist)
                    max_dist = dist;
            }

            for (i = 0; i < num_channels; i++)
            {
                vis_rgb_channels[i](y, x) = max_dist*255.0;
            }

            identity_score(y, x) = identity_score(y, x)*max_dist;
        }
    }

    errno_t err = write_png_file(vis_image_path.c_str(), vis_rgb_channels, 3);
}

void calculate_average_intensity(MatrixXd& intensity, int x0, int y0, int wx, int wy, double& average_intensity)
{
    int x, y;

    int width = intensity.cols();
    int height = intensity.rows();
    average_intensity = 0.0;
    int wsize = (2 * wx + 1)*(2 * wy + 1);

    if (y0 < wy || y0 + wy >= height || x0 < wx || x0 + wx >= width)
        return;

    for (y = y0 - wy; y <= y0 + wy; y++) {
        for (x = x0 - wx; x <= x0 + wx; x++) {
            average_intensity += intensity(y, x);
        }
    }

    average_intensity /= wsize;
}

void compare_average_intensities(double average_intensity_1, std::vector<std::vector<double>>& average_intensity_2, int num_channels, std::string vis_image_path, MatrixXd& identity_score)
{
    int i;
    int x, y;
    double dist;
    MatrixXd vis_rgb_channels[3];
    int height = average_intensity_2.size();
    int width = average_intensity_2[0].size();

    for (i = 0; i < num_channels; i++)
    {
        vis_rgb_channels[i] = MatrixXd::Zero(height, width);
    }

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            dist = abs(average_intensity_1 - average_intensity_2[y][x]);

            for (i = 0; i < num_channels; i++)
            {
                vis_rgb_channels[i](y, x) = dist;
            }

            identity_score(y, x) = identity_score(y, x) * dist / 255.0;
        }
    }

    errno_t err = write_png_file(vis_image_path.c_str(), vis_rgb_channels, 3);
}


void visualize_identity_score(MatrixXd& identity_score, int num_channels, std::string vis_image_path)
{
    int i;
    int x, y;
    double dist;
    MatrixXd vis_rgb_channels[3];

    int width = identity_score.cols();
    int height = identity_score.rows();


    for (i = 0; i < num_channels; i++)
    {
        vis_rgb_channels[i] = MatrixXd::Zero(height, width);
    }

    int min_x, min_y;
    identity_score.minCoeff(&min_y, &min_x);

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            if (x == min_x && y == min_y) 
            {
                vis_rgb_channels[0](y, x) = 255.0;
            }
            else
            {
                for (i = 0; i < num_channels; i++)
                {
                    vis_rgb_channels[i](y, x) = identity_score(y, x) * 255.0;
                }
            }
        }
    }

    errno_t err = write_png_file(vis_image_path.c_str(), vis_rgb_channels, 3);
}

void calculate_edge_score(MatrixXd& rgb_channel, int x0, int y0, double& edge_score)
{
    int bin;
    int x, y;
    int wx = 1;
    int wy = 1;

    edge_score = 0.0;
    int width = rgb_channel.cols();
    int height = rgb_channel.rows();

    if (y0 < wy || y0 + wy >= height || x0 < wx || x0 + wx >= width)
        return;

    x = x0;
    y = y0;

    double dfdx = 
        (rgb_channel(y - 1, x + 1) - rgb_channel(y - 1, x - 1) +
        rgb_channel(y, x + 1) - rgb_channel(y, x - 1) +
        rgb_channel(y + 1, x + 1) - rgb_channel(y + 1, x - 1)) / 6.0;

    double dfdy =
        (rgb_channel(y + 1, x - 1) - rgb_channel(y - 1, x - 1) +
            rgb_channel(y + 1, x) - rgb_channel(y - 1, x) +
            rgb_channel(y + 1, x + 1) - rgb_channel(y - 1, x + 1)) / 6.0;

    edge_score = sqrt(dfdx*dfdx + dfdy*dfdy);
}

void compare_edge_scores(double edge_score_1[], std::vector<std::vector<double>> edge_score_2[], int num_channels, std::string vis_image_path, MatrixXd& identity_score)
{
    int i;
    int x, y;
    double diff, max_diff;
    MatrixXd vis_rgb_channels[3];
    MatrixXd edge_diff;
    int height = edge_score_2[0].size();
    int width = edge_score_2[0][0].size();

    for (i = 0; i < num_channels; i++)
    {
        vis_rgb_channels[i] = MatrixXd::Zero(height, width);
        edge_diff = MatrixXd::Zero(height, width);
    }

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            max_diff = 0.0;
            for (i = 0; i < num_channels; i++)
            {

                diff = abs(edge_score_1[i] - edge_score_2[i][y][x]);

                if (max_diff < diff)
                    max_diff = diff;
            }

            edge_diff(y, x) = max_diff;
        }
    }


    double maxCoeff;
    maxCoeff = edge_diff.maxCoeff();
    edge_diff /= maxCoeff;

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            for (i = 0; i < num_channels; i++)
            {
                vis_rgb_channels[i](y, x) = edge_diff(y,x)*255.0;
            }

            identity_score(y, x) = identity_score(y, x)*edge_diff(y, x);
        }
    }

    errno_t err = write_png_file(vis_image_path.c_str(), vis_rgb_channels, 3);
}

int main(int argc, char **argv)
{
    MatrixXd rgb_channels_1[3], rgb_channels_2[3], rgb_channels_3[3];
    MatrixXd intensity_1, intensity_2;
    MatrixXd identity_score;
    int num_channels = 3;

    std::string image_path[2];
    std::string identity_score_image_path;
    std::string edge_image_path;
    std::string average_intensity_image_path;
    std::string histogram_diff_image_path;
    std::string arg;
    int x0, y0, wx, wy;
    int num_bins;
    int x, y;

/*
    Command line arguments:
-num_bins 10 -wx 3 -wy 3 -x0 75 -y0 37 -i1 C:\Users\Sam\Documents\EyeSignalsProjects\images\image1_025.png -i2 C:\Users\Sam\Documents\EyeSignalsProjects\images\image2_025.png -ohist C:\Users\Sam\Documents\EyeSignalsProjects\images\histogram_diff.png -oavint C:\Users\Sam\Documents\EyeSignalsProjects\images\av_intensity_diff.png -oid C:\Users\Sam\Documents\EyeSignalsProjects\images\identity_score.png -oedge C:\Users\Sam\Documents\EyeSignalsProjects\images\edge_diff.png
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
        if (arg == "-oid")
        {
            i++;
            identity_score_image_path = argv[i];
        }
        if (arg == "-oedge")
        {
            i++;
            edge_image_path = argv[i];
        }

        if (arg == "-ohist")
        {
            i++;
            histogram_diff_image_path = argv[i];
        }

        if (arg == "-oavint")
        {
            i++;
            average_intensity_image_path = argv[i];
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

    preprocess_channels(rgb_channels_1, num_channels, intensity_1);
    preprocess_channels(rgb_channels_2, num_channels, intensity_2);

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

    int width = rgb_channels_2[0].cols();
    int height = rgb_channels_2[0].rows();

    identity_score = MatrixXd::Ones(height, width);


    // histogram statistic
    VectorXd histogram_1[3];
    std::vector<std::vector<VectorXd>> histograms_2[3];

    for (i = 0; i < 3; i++)
    {
        calculate_histogram(rgb_channels_1[i], x0, y0, wx, wy, histogram_1[i], num_bins);
    }

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

    compare_histograms(histogram_1, histograms_2, num_channels, histogram_diff_image_path, identity_score);

    // on-edge statistic
    double edge_score_1[3];
    std::vector<std::vector<double>> edge_score_2[3];

    for (i = 0; i < 3; i++)
    {
        calculate_edge_score(rgb_channels_1[i], x0, y0, edge_score_1[i]);
    }

    for (i = 0; i < 3; i++)
    {
        edge_score_2[i].resize(height);
        for (y = 0; y < height; y++) {
            edge_score_2[i][y].resize(width);
            for (x = 0; x < width; x++) {
                calculate_edge_score(rgb_channels_2[i], x, y, edge_score_2[i][y][x]);
            }
        }
    }

    compare_edge_scores(edge_score_1, edge_score_2, num_channels, edge_image_path, identity_score);

#ifdef use_average_intensity
    // average intensity statistic
    double average_intensity_1;
    std::vector<std::vector<double>> average_intensity_2;

    calculate_average_intensity(intensity_1, x0, y0, wx, wy, average_intensity_1);
    average_intensity_2.resize(height);
    for (y = 0; y < height; y++) {
        average_intensity_2[y].resize(width);
        for (x = 0; x < width; x++) {
            calculate_average_intensity(intensity_2, x, y, wx, wy, average_intensity_2[y][x]);
        }
    }
    compare_average_intensities(average_intensity_1, average_intensity_2, num_channels, average_intensity_image_path, identity_score);
#endif

    visualize_identity_score(identity_score, num_channels, identity_score_image_path);

    return 0;
}

