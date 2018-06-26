#pragma once

class param_context
{
public:
    param_context()
    {
        x1 = 0;
        y1 = 0;
        batch_mode = false;
        num_candidates = 3;
        num_bins = 10;
        x0 = 0;
        y0 = 0;
        num_channels = 3;
        cand_coef = 0.1;
        draw_corr_point = true;
        wx = 3;
        wy = 3;
    }

    /*
    Command line arguments:
    -x0 136 -y0 53 -x1 90 -y1 54 -num_bins 10 -wx 3 -wy 3 -i1 C:\Users\Sam\Documents\EyeSignalsProjects\images\image1_025.png -i2 C:\Users\Sam\Documents\EyeSignalsProjects\images\image2_025.png -oid C:\Users\Sam\Documents\EyeSignalsProjects\images\identity_score.png -oavint C:\Users\Sam\Documents\EyeSignalsProjects\images\av_intensity_diff.png -ohist C:\Users\Sam\Documents\EyeSignalsProjects\images\histogram_diff.png -ograd C:\Users\Sam\Documents\EyeSignalsProjects\images\gradient_diff.png -wdif C:\Users\Sam\Documents\EyeSignalsProjects\images\win_diff.png
    -batch  -candcoef 0.01 -numcand 3 -num_bins 10 -wx 5 -wy 5 -i1 C:\Users\Sam\Documents\EyeSignalsProjects\images\image1_025.png -i2 C:\Users\Sam\Documents\EyeSignalsProjects\images\image2_025.png
    -x0 102 -y0 52 -x1 58 -y1 51 -candcoef 0.01 -num_bins 10 -wx 3 -wy 3 -i1 C:\Users\Sam\Documents\EyeSignalsProjects\images\image1_025.png -i2 C:\Users\Sam\Documents\EyeSignalsProjects\images\image2_025.png -oid C:\Users\Sam\Documents\EyeSignalsProjects\images\identity_score.png -oavint C:\Users\Sam\Documents\EyeSignalsProjects\images\av_intensity_diff.png -ohist C:\Users\Sam\Documents\EyeSignalsProjects\images\histogram_diff.png -ograd C:\Users\Sam\Documents\EyeSignalsProjects\images\gradient_diff.png -ovline C:\Users\Sam\Documents\EyeSignalsProjects\images\vline_diff.png
    -x0 78 -y0 12 -x1 32 -y1 8 -num_bins 10 -wx 3 -wy 3 -i1 C:\Users\Sam\Documents\EyeSignalsProjects\images\image1_025.png -i2 C:\Users\Sam\Documents\EyeSignalsProjects\images\image2_025.png -oid C:\Users\Sam\Documents\EyeSignalsProjects\images\identity_score.png -oavint C:\Users\Sam\Documents\EyeSignalsProjects\images\av_intensity_diff.png -ohist C:\Users\Sam\Documents\EyeSignalsProjects\images\histogram_diff.png -ograd C:\Users\Sam\Documents\EyeSignalsProjects\images\gradient_diff.png
    */

    void get_params_from_command_line(int argc, char **argv)
    {
        std::string arg;

        int i = 1;
        while (i < argc)
        {
            arg = argv[i];

            if (arg == "-vis_grad")
            {
                i++;
                vis_gradient_image_path = argv[i];
            }

            if (arg == "-candcoef")
            {
                i++;
                cand_coef = atof(argv[i]);
            }

            if (arg == "-nocorpnt")
            {
                draw_corr_point = false;
            }

            if (arg == "-batch")
            {
                batch_mode = true;
            }
            if (arg == "-numcand")
            {
                i++;
                num_candidates = atoi(argv[i]);
            }
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
            if (arg == "-ograd")
            {
                i++;
                gradient_diff_image_path = argv[i];
            }

            if (arg == "-ohist")
            {
                i++;
                histogram_diff_image_path = argv[i];
            }

            if (arg == "-wdif")
            {
                i++;
                windows_diff_image_path = argv[i];
            }

            if (arg == "-ovline")
            {
                i++;
                vline_diff_image_path = argv[i];
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

            if (arg == "-x1")
            {
                i++;
                x1 = atoi(argv[i]);
            }

            if (arg == "-y1")
            {
                i++;
                y1 = atoi(argv[i]);
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
    }

    bool batch_mode;
    int num_candidates;
    std::string image_path[2];
    std::string identity_score_image_path;
    std::string gradient_diff_image_path;
    std::string vis_gradient_image_path;
    std::string average_intensity_image_path;
    std::string histogram_diff_image_path;
    std::string vline_diff_image_path;
    std::string windows_diff_image_path;
    
    int x0;  // a point of the image1 for which correspondence is calculated in the image2
    int y0;
    int x1;  // a point of the image2 which is the ground truth correspondence for (x0,y0)
    int y1;
    int wx;  // half window size in which a feature is calculated
    int wy;
    int num_bins;
    int num_channels;
    double cand_coef;
    bool draw_corr_point;
};
