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
        max_diff = 0.1;
        num_chan_ranges = 10;
        statistic_localization = 0.1;
        learn_statistic_iterations = 20;
        num_ransack_iterations = 20;
        F_err_thresh = 2*2*4;
        statistic_type = 0;
        num_displayed_correspondences = 20;
        pass_ratio_thresh = 0.1;
        seed_random_engines = false;
        num_stat_components = 3;
    }

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

            if (arg == "-ferr")
            {
                i++;
                F_err_thresh = atof(argv[i]);
            }

            if (arg == "-stloc")
            {
                i++;
                statistic_localization = atof(argv[i]);
            }
            
            if (arg == "-nocorpnt")
            {
                draw_corr_point = false;
            }

            if (arg == "-ndc")
            {
                i++;
                num_displayed_correspondences = atoi(argv[i]);
            }

            if (arg == "-seed")
            {
                seed_random_engines = true;
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
            if (arg == "-nri")
            {
                i++;
                num_ransack_iterations = atoi(argv[i]);
            }
            if (arg == "-statt")
            {
                i++;
                statistic_type = atoi(argv[i]);
            }

            if (arg == "-prt")
            {
                i++;
                pass_ratio_thresh = atof(argv[i]);
            }
            
            if (arg == "-nsc")
            {
                i++;
                num_stat_components = atoi(argv[i]);
            }

            if (arg == "-ncr")
            {
                i++;
                num_chan_ranges = atoi(argv[i]);
            }

            if (arg == "-lsi")
            {
                i++;
                learn_statistic_iterations = atoi(argv[i]);
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
            if (arg == "-ocor")
            {
                i++;
                correspondence_image_path = argv[i];
            }

            if (arg == "-opcor")
            {
                i++;
                point_correspondence_image_path = argv[i];
            }

            if (arg == "-ogdst")
            {
                i++;
                good_statistics_file_path = argv[i];
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

            if (arg == "-max_diff")
            {
                i++;
                max_diff = atof(argv[i]);
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
    std::string correspondence_image_path;
    std::string good_statistics_file_path;
    std::string point_correspondence_image_path;
    
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
    double max_diff;
    int num_chan_ranges;
    double statistic_localization;
    int learn_statistic_iterations;
    int num_ransack_iterations;
    double  F_err_thresh;
    int statistic_type;
    int num_displayed_correspondences;
    double pass_ratio_thresh;
    bool seed_random_engines;
    int num_stat_components;
};
