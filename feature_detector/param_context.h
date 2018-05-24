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
    }

    /*
    Command line arguments:
    -batch -numcand 3 -num_bins 10 -wx 3 -wy 3 -i1 C:\Users\Sam\Documents\EyeSignalsProjects\images\image1_025.png -i2 C:\Users\Sam\Documents\EyeSignalsProjects\images\image2_025.png
    -x0 136 -y0 53 -x1 90 -y1 54 -num_bins 10 -wx 3 -wy 3 -i1 C:\Users\Sam\Documents\EyeSignalsProjects\images\image1_025.png -i2 C:\Users\Sam\Documents\EyeSignalsProjects\images\image2_025.png -oid C:\Users\Sam\Documents\EyeSignalsProjects\images\identity_score.png -oavint C:\Users\Sam\Documents\EyeSignalsProjects\images\av_intensity_diff.png -ohist C:\Users\Sam\Documents\EyeSignalsProjects\images\histogram_diff.png -ograd C:\Users\Sam\Documents\EyeSignalsProjects\images\gradient_diff.png
    -x0 78 -y0 12 -x1 32 -y1 8 -num_bins 10 -wx 3 -wy 3 -i1 C:\Users\Sam\Documents\EyeSignalsProjects\images\image1_025.png -i2 C:\Users\Sam\Documents\EyeSignalsProjects\images\image2_025.png -oid C:\Users\Sam\Documents\EyeSignalsProjects\images\identity_score.png -oavint C:\Users\Sam\Documents\EyeSignalsProjects\images\av_intensity_diff.png -ohist C:\Users\Sam\Documents\EyeSignalsProjects\images\histogram_diff.png -ograd C:\Users\Sam\Documents\EyeSignalsProjects\images\gradient_diff.png
    */

    void get_params_from_command_line(int argc, char **argv)
    {
        std::string arg;

        int i = 1;
        while (i < argc)
        {
            arg = argv[i];
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
                gradient_image_path = argv[i];
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

    int x1;
    int y1;
    bool batch_mode;
    int num_candidates;
    std::string image_path[2];
    std::string identity_score_image_path;
    std::string gradient_image_path;
    std::string average_intensity_image_path;
    std::string histogram_diff_image_path;
    int x0;
    int y0;
    int wx;
    int wy;
    int num_bins;
    int num_channels;
};
