// feature_detector.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "lpngwrapper.hpp"
#include "param_context.h"
#include <list>
#include <vector>
#include <iostream>

#include "image_correspondence.h"

// version 5-24-18-2-03

// coordinates of the corresponidng point on the second image
void usage(void)
{
    fprintf(stderr, "usage: feature_detector <input_image_1_file_path> <input_image_2_file_path> <output_image_3_file_path>\n");
    exit(1);
}

// input: in [0,1]
// output: in [0.5,1]  
double dist_func(double in_val)
{
    double in_min_bound = 0.0;
    double in_max_bound = 1.0;
    double out_min_bound = 0.8;
    double out_max_bound = 1.0;
    double out_val = (in_val - in_min_bound) / (in_max_bound - in_min_bound);
    out_val = out_min_bound + out_val*(out_max_bound - out_min_bound);
    if (out_val < out_min_bound)
        out_val = out_min_bound;
    if (out_val > out_max_bound)
        out_val = out_max_bound;
    return out_val;
}

void preprocess_channels(param_context& pc, MatrixXd rgb_channels[], MatrixXd& intensity)
{
    int x, y;
    double r, g, b, sum;

    int width = rgb_channels[0].cols();
    int height = rgb_channels[0].rows();
    intensity = MatrixXd::Zero(height, width);

    assert(pc.num_channels == 3);

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
            intensity(y, x) = sum / 3.0 / 255.0;
        }
    }
}

void DrawCorrespondingPoint(param_context& pc, MatrixXd vis_rgb_channels[])
{
    int i;

    if (!pc.draw_corr_point)
        return;

    for (i = 0; i < pc.num_channels; i++)
    {
        vis_rgb_channels[i](pc.y1, pc.x1) = 0;
    }
    vis_rgb_channels[0](pc.y1, pc.x1) = 255.0;
}

void calculate_histogram(param_context& pc, MatrixXd& rgb_channel, int x0, int y0, VectorXd& hist)
{
    int bin;
    int x, y;

    hist = VectorXd::Zero(pc.num_bins);
    int width = rgb_channel.cols();
    int height = rgb_channel.rows();
    double bin_size = 1.0 / pc.num_bins;

    if (y0 < pc.wy || y0 + pc.wy >= height || x0 < pc.wx || x0 + pc.wx >= width)
        return;

    for (y = y0 - pc.wy; y <= y0 + pc.wy; y++) {
        for (x = x0 - pc.wx; x <= x0 + pc.wx; x++) {
            bin = rgb_channel(y, x) / bin_size;

            if (bin >= pc.num_bins)
                bin = pc.num_bins - 1;
            hist(bin)++;
        }
    }

    hist.normalize();

}

void compare_histograms(param_context& pc, VectorXd histogram_1[], std::vector<std::vector<VectorXd>> histograms_2[], MatrixXd& identity_score)
{
    int i;
    int x, y;
    double dist, max_dist;
    MatrixXd vis_rgb_channels[3];
    int height = histograms_2[0].size();
    int width = histograms_2[0][0].size();

    for (i = 0; i < pc.num_channels; i++)
    {
        vis_rgb_channels[i] = MatrixXd::Zero(height, width);
    }

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            max_dist = 0.0;
            for (i = 0; i < pc.num_channels; i++)
            {
                dist = (histogram_1[i] - histograms_2[i][y][x]).cwiseAbs().maxCoeff();
                if (max_dist < dist)
                    max_dist = dist;
            }

            for (i = 0; i < pc.num_channels; i++)
            {
                vis_rgb_channels[i](y, x) = max_dist*255.0;
            }

            identity_score(y, x) = identity_score(y, x)*dist_func(max_dist);
        }
    }

    if (pc.batch_mode)
        return;

    DrawCorrespondingPoint(pc, vis_rgb_channels);
    errno_t err = write_png_file(pc.histogram_diff_image_path.c_str(), vis_rgb_channels, 3);
}

void calculate_average_intensity(param_context& pc, MatrixXd& intensity, int x0, int y0, int wx, int wy, double& average_intensity)
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

void compare_average_intensities(param_context& pc, double average_intensity_1, std::vector<std::vector<double>>& average_intensity_2, MatrixXd& identity_score)
{
    int i;
    int x, y;
    double dist;
    double max_intensity;
    MatrixXd vis_rgb_channels[3];
    int height = average_intensity_2.size();
    int width = average_intensity_2[0].size();

    for (i = 0; i < pc.num_channels; i++)
    {
        vis_rgb_channels[i] = MatrixXd::Zero(height, width);
    }

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            max_intensity = std::max(average_intensity_1, average_intensity_2[y][x]);
#define intensity_threshold 0.01
            if (max_intensity < intensity_threshold)
                dist = 0;
            else
                dist = abs(average_intensity_1 - average_intensity_2[y][x]) / max_intensity;

            for (i = 0; i < pc.num_channels; i++)
            {
                vis_rgb_channels[i](y, x) = dist*255.0;
            }

            identity_score(y, x) = identity_score(y, x) * dist_func(dist);
        }
    }

    if (pc.batch_mode)
        return;

    DrawCorrespondingPoint(pc, vis_rgb_channels);
    errno_t err = write_png_file(pc.average_intensity_image_path.c_str(), vis_rgb_channels, 3);
}


void visualize_identity_score(param_context& pc, MatrixXd& identity_score, int xb, int yb)
{
    int i;
    int x, y;
    double dist;
    MatrixXd vis_rgb_channels[3];

    int width = identity_score.cols();
    int height = identity_score.rows();


    for (i = 0; i < pc.num_channels; i++)
    {
        vis_rgb_channels[i] = MatrixXd::Zero(height, width);
    }

    double minCoeff, maxCoeff;
    minCoeff = identity_score.minCoeff();
    maxCoeff = identity_score.maxCoeff();
    double candidateCoeff = minCoeff + (maxCoeff - minCoeff)*pc.cand_coef;
    int candidateCoeffCnt = 0;
    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            double score = identity_score(y, x);
            if (!pc.batch_mode)
            {
                if (score < candidateCoeff)
                {
                    if (x == pc.x1 && y == pc.y1)
                        vis_rgb_channels[1](y, x) = 255.0;
                    else
                        vis_rgb_channels[0](y, x) = 255.0;
                }
                else
                {
                    if (x == pc.x1 && y == pc.y1)
                        vis_rgb_channels[2](y, x) = 255.0;
                    else
                    {
                        for (i = 0; i < pc.num_channels; i++)
                        {
                            vis_rgb_channels[i](y, x) = identity_score(y, x) * 255.0;
                        }
                    }
                }
            }
            else
            {
                if (score < candidateCoeff)
                    candidateCoeffCnt++;
            }
        }
    }

    if (pc.batch_mode)
    {
        if (candidateCoeffCnt <= pc.num_candidates)
            std::cout << "candidateCoeffCnt=" << candidateCoeffCnt << " xb=" << xb << " yb=" << yb << "\n";
        return;
    }

    errno_t err = write_png_file(pc.identity_score_image_path.c_str(), vis_rgb_channels, 3);
}

void calculate_gradient(MatrixXd& rgb_channel, int x0, int y0, double& gradient)
{
    int bin;
    int x, y;
    int wx = 1;
    int wy = 1;

    gradient = 0.0;
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

    gradient = sqrt(dfdx*dfdx + dfdy*dfdy);
}

void compare_gradients(param_context& pc, double gradient_1[], std::vector<std::vector<double>> gradient_2[], MatrixXd& identity_score)
{
    int i;
    int x, y;
    double diff, max_diff;
    double max_gradient;
    double maxCoeff;

    MatrixXd vis_rgb_channels[3];
    MatrixXd gradient_diff;
    int height = gradient_2[0].size();
    int width = gradient_2[0][0].size();

    assert(pc.num_channels == 3);

    if (pc.vis_gradient_image_path.length() > 0)
    {
        MatrixXd vis_gradient[3];
        MatrixXd max_gradient = MatrixXd::Zero(height, width);

        for (y = 0; y < height; y++)
        {
            for (x = 0; x < width; x++)
            {
                double max_val = 0.0;
                for (i = 0; i < pc.num_channels; i++)
                {
                    if (max_val < gradient_2[i][y][x])
                        max_val = gradient_2[i][y][x];
                }

                max_gradient(y,x) = max_val;
            }
        }
        maxCoeff = max_gradient.maxCoeff();
        max_gradient /= maxCoeff;
        max_gradient *= 255.0;
        for (i = 0; i < pc.num_channels; i++)
        {
            vis_gradient[i] = max_gradient;
        }
        errno_t err = write_png_file(pc.vis_gradient_image_path.c_str(), vis_gradient, 3);
    }

    for (i = 0; i < pc.num_channels; i++)
    {
        vis_rgb_channels[i] = MatrixXd::Zero(height, width);
        gradient_diff = MatrixXd::Zero(height, width);
    }

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            max_diff = 0.0;
            for (i = 0; i < pc.num_channels; i++)
            {
                max_gradient = std::max(gradient_1[i], gradient_2[i][y][x]);
#define gradient_threshold 0.001
                if (max_gradient > gradient_threshold)
                {
                    diff = abs(gradient_1[i] - gradient_2[i][y][x]) / max_gradient;
                    if (max_diff < diff)
                        max_diff = diff;
                }
            }

            gradient_diff(y, x) = max_diff;
        }
    }

    maxCoeff = gradient_diff.maxCoeff();
    gradient_diff /= maxCoeff;

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            for (i = 0; i < pc.num_channels; i++)
            {
                vis_rgb_channels[i](y, x) = gradient_diff(y, x)*255.0;
            }

            identity_score(y, x) = identity_score(y, x)*dist_func(gradient_diff(y, x));
        }
    }

    if (pc.batch_mode)
        return;

    DrawCorrespondingPoint(pc, vis_rgb_channels);
    errno_t err = write_png_file(pc.gradient_diff_image_path.c_str(), vis_rgb_channels, 3);
}

void calculate_vline(param_context& pc, MatrixXd& rgb_channel, int x0, int y0, VectorXd& vline)
{
    int b;
    int x, y;

    vline = VectorXd::Zero(2*pc.wy+1);
    int width = rgb_channel.cols();
    int height = rgb_channel.rows();

    if (y0 < pc.wy || y0 + pc.wy >= height)
        return;

    x = x0;
    b = 0;
    for (y = y0 - pc.wy; y <= y0 + pc.wy; y++) {
        vline(b) = rgb_channel(y, x);
        b++;
    }
}

void compare_vlines(param_context& pc, VectorXd vline_1[], std::vector<std::vector<VectorXd>> vlines_2[], MatrixXd& identity_score)
{
    int i;
    int x, y;
    double dist, max_dist;
    MatrixXd vis_rgb_channels[3];
    int height = vlines_2[0].size();
    int width = vlines_2[0][0].size();

    for (i = 0; i < pc.num_channels; i++)
    {
        vis_rgb_channels[i] = MatrixXd::Zero(height, width);
    }

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            max_dist = 0.0;
            for (i = 0; i < pc.num_channels; i++)
            {
                dist = (vline_1[i] - vlines_2[i][y][x]).cwiseAbs().maxCoeff();
                if (max_dist < dist)
                    max_dist = dist;
            }

            for (i = 0; i < pc.num_channels; i++)
            {
                vis_rgb_channels[i](y, x) = max_dist*255.0;
            }

            identity_score(y, x) = identity_score(y, x)*dist_func(max_dist);
        }
    }

    if (pc.batch_mode)
        return;

    DrawCorrespondingPoint(pc, vis_rgb_channels);
    errno_t err = write_png_file(pc.vline_diff_image_path.c_str(), vis_rgb_channels, 3);
}

bool compare_windows(param_context& pc, MatrixXd rgb_channels_1[], MatrixXd& intensity_1, MatrixXd rgb_channels_2[], MatrixXd& intensity_2, int x0, int y0, int x1, int y1)
{
    int x, y, c;

    double avg, diff;
    double max_diff = pc.max_diff;   // 0.15; // 0.3;  // 0.28; 
    int width = rgb_channels_1[0].cols();
    int height = rgb_channels_1[0].rows();

    if (y0 < pc.wy || y0 + pc.wy >= height || x0 < pc.wx || x0 + pc.wx >= width)
        return false;
    if (y1 < pc.wy || y1 + pc.wy >= height || x1 < pc.wx || x1 + pc.wx >= width)
        return false;

    for (y = - pc.wy; y <= pc.wy; y++) {
        for (x = - pc.wx; x <= pc.wx; x++) {
                for (c = 0; c < pc.num_channels; c++) {
                    avg = (rgb_channels_1[c](y + y0, x + x0) + rgb_channels_2[c](y + y1, x + x1))/2;
                    diff = rgb_channels_1[c](y + y0, x + x0) - rgb_channels_2[c](y + y1, x + x1);
                    if (avg > 0 && diff / avg > max_diff)
                    {
                        return false;
                    }
                }

                avg = (intensity_1(y + y0, x + x0) + intensity_2(y + y1, x + x1))/2;
                diff = abs(intensity_1(y + y0, x + x0) - intensity_2(y + y1, x + x1));
                if (avg > 0 && diff / avg > max_diff)
                {
                    return false;
                }
        }
    }

    return true;
}

void compare_images(param_context& pc, MatrixXd rgb_channels_1[], MatrixXd& intensity_1,  MatrixXd rgb_channels_2[], MatrixXd& intensity_2, MatrixXd& identity_score)
{
    int i;
    int x, y;
    MatrixXd vis_rgb_channels[3];
    int width = rgb_channels_1[0].cols();
    int height = rgb_channels_1[0].rows();
    MatrixXd region;

    for (i = 0; i < pc.num_channels; i++)
    {
        vis_rgb_channels[i] = MatrixXd::Zero(height, width);
    }

    for (y = 0; y < height; y++)
    {
        for (x = 0; x < width; x++)
        {
            if (compare_windows(pc, rgb_channels_1, intensity_1, rgb_channels_2, intensity_2, pc.x0, pc.y0, x, y))
            {
                for (i = 0; i < pc.num_channels; i++)
                {
                    vis_rgb_channels[i](y, x) = 255.0;
                }
            }
            else
            {
                identity_score(y, x) = 0;
            }
        }
    }

    if (pc.batch_mode)
        return;

    DrawCorrespondingPoint(pc, vis_rgb_channels);
    errno_t err = write_png_file(pc.windows_diff_image_path.c_str(), vis_rgb_channels, 3);
}

/*
example command parameters:
-mnc 32 -prt 0.3 -nsc 2 -stloc 0.1 -ndc 40 -ncr 10 -statt 0 -ferr 16 -lsi 100 -nri 200 -max_diff 0.5 -x0 276 -y0 50 -x1 0 -y1 0 -num_bins 10 -wx 7 -wy 7 -i1 C:\Users\Sam\Documents\EyeSignalsProjects\images\image1.png -i2 C:\Users\Sam\Documents\EyeSignalsProjects\images\image2.png -oid C:\Users\Sam\Documents\EyeSignalsProjects\images\identity_score.png -oavint C:\Users\Sam\Documents\EyeSignalsProjects\images\av_intensity_diff.png -ohist C:\Users\Sam\Documents\EyeSignalsProjects\images\histogram_diff.png -ograd C:\Users\Sam\Documents\EyeSignalsProjects\images\gradient_diff.png -wdif C:\Users\Sam\Documents\EyeSignalsProjects\images\win_diff.png -ocor  C:\Users\Sam\Documents\EyeSignalsProjects\images\correspondence.png -ogdst  C:\Users\Sam\Documents\EyeSignalsProjects\images\good_statistics.txt  -opcor  C:\Users\Sam\Documents\EyeSignalsProjects\images\point_correspondence.png
-nsc 3 -seed -stloc 0.1 -prt 0.1 -ndc 20 -ncr 10 -statt 0 -ferr 16 -lsi 10 -nri 10 -max_diff 0.5 -x0 0 -y0 0 -x1 0 -y1 0 -num_bins 10 -wx 3 -wy 3 -i1 C : \Users\Sam\Documents\EyeSignalsProjects\images\image1.png -i2 C : \Users\Sam\Documents\EyeSignalsProjects\images\image2.png -oid C : \Users\Sam\Documents\EyeSignalsProjects\images\identity_score.png -oavint C : \Users\Sam\Documents\EyeSignalsProjects\images\av_intensity_diff.png -ohist C : \Users\Sam\Documents\EyeSignalsProjects\images\histogram_diff.png -ograd C : \Users\Sam\Documents\EyeSignalsProjects\images\gradient_diff.png -wdif C : \Users\Sam\Documents\EyeSignalsProjects\images\win_diff.png -ocor  C : \Users\Sam\Documents\EyeSignalsProjects\images\correspondence.png -ogdst  C : \Users\Sam\Documents\EyeSignalsProjects\images\good_statistics.txt -opcor  C : \Users\Sam\Documents\EyeSignalsProjects\images\point_correspondence.png

-nsc 1 -stloc 0.1 -ndc 20 -ncr 10 -statt 0 -ferr 16 -lsi 10 -nri 10 -max_diff 0.5 -x0 276 -y0 50 -x1 0 -y1 0 -num_bins 10 -wx 5 -wy 5 -i1 C:\Users\Sam\Documents\EyeSignalsProjects\images\image1_025.png -i2 C:\Users\Sam\Documents\EyeSignalsProjects\images\image2_025.png -oid C:\Users\Sam\Documents\EyeSignalsProjects\images\identity_score.png -oavint C:\Users\Sam\Documents\EyeSignalsProjects\images\av_intensity_diff.png -ohist C:\Users\Sam\Documents\EyeSignalsProjects\images\histogram_diff.png -ograd C:\Users\Sam\Documents\EyeSignalsProjects\images\gradient_diff.png -wdif C:\Users\Sam\Documents\EyeSignalsProjects\images\win_diff.png -ocor  C:\Users\Sam\Documents\EyeSignalsProjects\images\correspondence.png -ogdst  C:\Users\Sam\Documents\EyeSignalsProjects\images\good_statistics.txt  -opcor  C:\Users\Sam\Documents\EyeSignalsProjects\images\point_correspondence.png

-prt 0.1 -nsc 1 -seed -stloc 0.1 -ndc 20 -ncr 20 -statt 0 -ferr 4 -lsi 10 -nri 100 -max_diff 0.5 -x0 276 -y0 50 -x1 0 -y1 0 -num_bins 10 -wx 5 -wy 5 -i1 C:\Users\Sam\Documents\EyeSignalsProjects\images\image1.png -i2 C:\Users\Sam\Documents\EyeSignalsProjects\images\image2.png -oid C:\Users\Sam\Documents\EyeSignalsProjects\images\identity_score.png -oavint C:\Users\Sam\Documents\EyeSignalsProjects\images\av_intensity_diff.png -ohist C:\Users\Sam\Documents\EyeSignalsProjects\images\histogram_diff.png -ograd C:\Users\Sam\Documents\EyeSignalsProjects\images\gradient_diff.png -wdif C:\Users\Sam\Documents\EyeSignalsProjects\images\win_diff.png -ocor  C:\Users\Sam\Documents\EyeSignalsProjects\images\correspondence.png -ogdst  C:\Users\Sam\Documents\EyeSignalsProjects\images\good_statistics.txt  -opcor  C:\Users\Sam\Documents\EyeSignalsProjects\images\point_correspondence.png

-prt 0.1 -nsc 1 -seed -stloc 0.1 -ndc 20 -ncr 20 -statt 0 -ferr 4 -lsi 10 -nri 100 -max_diff 0.5 -x0 276 -y0 50 -x1 0 -y1 0 -num_bins 10 -wx 5 -wy 5 -i1 C:\Users\Sam\Documents\EyeSignalsProjects\images\image1.png -i2 C:\Users\Sam\Documents\EyeSignalsProjects\images\image2.png -oid C:\Users\Sam\Documents\EyeSignalsProjects\images\identity_score.png -oavint C:\Users\Sam\Documents\EyeSignalsProjects\images\av_intensity_diff.png -ohist C:\Users\Sam\Documents\EyeSignalsProjects\images\histogram_diff.png -ograd C:\Users\Sam\Documents\EyeSignalsProjects\images\gradient_diff.png -wdif C:\Users\Sam\Documents\EyeSignalsProjects\images\win_diff.png -ocor  C:\Users\Sam\Documents\EyeSignalsProjects\images\correspondence.png -ogdst  C:\Users\Sam\Documents\EyeSignalsProjects\images\good_statistics.txt  -opcor  C:\Users\Sam\Documents\EyeSignalsProjects\images\point_correspondence.png
*/
int main(int argc, char **argv)
{
    param_context pc;
    int i;
    int x, y;

    MatrixXd rgb_channels_1[3], rgb_channels_2[3], rgb_channels_3[3];
    MatrixXd intensity_1, intensity_2;
    MatrixXd identity_score;

    pc.get_params_from_command_line(argc, argv);

    c_image_correspondence image_correspondence(pc);
    image_correspondence.learn_statistic_parameters();

    read_png_file(pc.image_path[0].c_str(), rgb_channels_1, pc.num_channels);
    read_png_file(pc.image_path[1].c_str(), rgb_channels_2, pc.num_channels);

    bool point_inside =
            (pc.y1 >= 0 && pc.y1 < rgb_channels_1[0].rows() &&
            pc.x1 >= 0 && pc.x1 < rgb_channels_1[0].cols()) &&
            (pc.y0 >= 0 && pc.y0 < rgb_channels_1[0].rows() &&
            pc.x0 >= 0 && pc.x0 < rgb_channels_1[0].cols());

    assert(point_inside);
    if (!point_inside)
        return 1;

    preprocess_channels(pc, rgb_channels_1, intensity_1);
    preprocess_channels(pc, rgb_channels_2, intensity_2);

    int width = rgb_channels_2[0].cols();
    int height = rgb_channels_2[0].rows();

    for (int yb = pc.wy; yb < height - pc.wy; yb++)
    {
        for (int xb = pc.wx; xb < width - pc.wx; xb++)
        {
            if (!pc.batch_mode)
            {
                xb = pc.x0;
                yb = pc.y0;
            }

            identity_score = MatrixXd::Ones(height, width);

            compare_images(pc, rgb_channels_1, intensity_1, rgb_channels_2, intensity_2, identity_score);
#if 0
            // vertical line statistic
            VectorXd vline_1[3];
            std::vector<std::vector<VectorXd>> vlines_2[3];

            for (i = 0; i < 3; i++)
            {
                calculate_vline(pc, rgb_channels_1[i], xb, yb, vline_1[i]);
            }

            for (i = 0; i < 3; i++)
            {
                vlines_2[i].resize(height);
                for (y = 0; y < height; y++) {
                    vlines_2[i][y].resize(width);
                    for (x = 0; x < width; x++) {
                        calculate_vline(pc, rgb_channels_2[i], x, y, vlines_2[i][y][x]);
                    }
                }
            }

            compare_vlines(pc, vline_1, vlines_2, identity_score);

            // histogram statistic
            VectorXd histogram_1[3];
            std::vector<std::vector<VectorXd>> histograms_2[3];

            for (i = 0; i < 3; i++)
            {
                calculate_histogram(pc, rgb_channels_1[i], xb, yb, histogram_1[i]);
            }

            for (i = 0; i < 3; i++)
            {
                histograms_2[i].resize(height);
                for (y = 0; y < height; y++) {
                    histograms_2[i][y].resize(width);
                    for (x = 0; x < width; x++) {
                        calculate_histogram(pc, rgb_channels_2[i], x, y, histograms_2[i][y][x]);
                    }
                }
            }

            compare_histograms(pc, histogram_1, histograms_2, identity_score);

            // gradient statistic
            double gradient_1[3];
            std::vector<std::vector<double>> gradient_2[3];

            for (i = 0; i < 3; i++)
            {
                calculate_gradient(rgb_channels_1[i], xb, yb, gradient_1[i]);
            }

            for (i = 0; i < 3; i++)
            {
                gradient_2[i].resize(height);
                for (y = 0; y < height; y++) {
                    gradient_2[i][y].resize(width);
                    for (x = 0; x < width; x++) {
                        calculate_gradient(rgb_channels_2[i], x, y, gradient_2[i][y][x]);
                    }
                }
            }

            compare_gradients(pc, gradient_1, gradient_2, identity_score);

            // average intensity statistic
            double average_intensity_1;
            std::vector<std::vector<double>> average_intensity_2;

            double av_int_wx = 1;
            double av_int_wy = 1;
            calculate_average_intensity(pc, intensity_1, xb, yb, av_int_wx, av_int_wy, average_intensity_1);
            average_intensity_2.resize(height);
            for (y = 0; y < height; y++) {
                average_intensity_2[y].resize(width);
                for (x = 0; x < width; x++) {
                    calculate_average_intensity(pc, intensity_2, x, y, av_int_wx, av_int_wy, average_intensity_2[y][x]);
                }
            }
            compare_average_intensities(pc, average_intensity_1, average_intensity_2, identity_score);
#endif
            visualize_identity_score(pc, identity_score, xb, yb);

            if (!pc.batch_mode)
                break;
        }
        if (!pc.batch_mode)
            break;
    }
    return 0;
}

