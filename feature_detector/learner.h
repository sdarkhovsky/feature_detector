#pragma once

#include "lpngwrapper.hpp"
#include "param_context.h"

#include <vector>
#include <map>
#include <random>
#include <iostream>
#include <fstream>  

using namespace std;

namespace Learner {

    class c_statistic
    {
    public:

        void calculate_linear_statistic_parameters(param_context& pc)
        {
            std::uniform_real_distribution<double> distribution(-1.0, 1.0);
            if (pc.seed_random_engines)
            {
                statistic_generator.seed(time(NULL));
            }

            int wx = pc.wx;
            int wy = pc.wy;

            int wsize = (2 * wx + 1)*(2 * wy + 1);

            params = MatrixXd::Zero(2 * wy + 1, 2 * wx + 1);
            for (int i1 = -wy; i1 <= wy; i1++)
            {
                for (int i2 = -wx; i2 <= wx; i2++)
                {
                    params(i1 + wy, i2 + wx) = distribution(statistic_generator);
                }
            }

            statistic_type = (c_statistic_type) pc.statistic_type;
            switch (statistic_type)
            {
            case c_statistic::differential:
            {
                // make total sum 0
                params(wy, wx) = -(params.sum() - params(wy, wx));
            }
            default:
            {
            }
            }
        }

        MatrixXd params;
        enum c_statistic_type { differential = 0, linear = 1, sigmoid = 2 };
        c_statistic_type statistic_type;

        std::default_random_engine statistic_generator;
    };

    class c_statistic_channel
    {
    public:
        void calculate_image_statistic(param_context& pc, const MatrixXd& image_channel, c_statistic& statistic)
        {
            int x, y, x0, y0;
            double sum;

            int width = image_channel.cols();
            int height = image_channel.rows();

            int wx = (statistic.params.cols() - 1) / 2;
            int wy = (statistic.params.rows() - 1) / 2;

            statistic_val = MatrixXd::Zero(height, width);

            for (y0 = wy; y0 < height - wy; y0++) {
                for (x0 = wx; x0 < width - wx; x0++) {
                    sum = 0.0;
                    for (y = y0 - wy; y <= y0 + wy; y++) {
                        for (x = x0 - wx; x <= x0 + wx; x++) {
                            sum += image_channel(y, x)*statistic.params(y - y0 + wy, x - x0 + wx);
                        }
                    }

                    switch (statistic.statistic_type)
                    {
                    case c_statistic::sigmoid:
                    {
                        statistic_val(y0, x0) = 1 / (1 + exp(-sum));
                    }
                    default:
                    {
                        statistic_val(y0, x0) = sum;
                    }
                    }
                }
            }

            minCoeff = statistic_val.minCoeff();
            maxCoeff = statistic_val.maxCoeff();
            range_length = (maxCoeff - minCoeff) / pc.num_chan_ranges;
        }

        MatrixXd statistic_val;
        double minCoeff, maxCoeff, range_length;
    };

    class c_range_key
    {
    public:
        bool operator() (const c_range_key& lhs, const c_range_key& rhs) const
        {
            int n = lhs.ranges.size();
            if (n > rhs.ranges.size())
                n = rhs.ranges.size();
            for (int i = 0; i < n; i++)
            {
                if (lhs.ranges[i] < rhs.ranges[i])
                    return true;
                if (lhs.ranges[i] > rhs.ranges[i])
                    return false;
            }
            return false;
        }

        vector<int> ranges;
    };

    class c_point
    {
    public:
        int x, y;
    };

    class c_point_set
    {
    public:
        c_point_set()
        {
        }

        void calculate_boundaries()
        {
            if (points.size() == 0)
                return;

            min_x = points[0].x;
            max_x = points[0].x;
            min_y = points[0].y;
            max_y = points[0].y;

            for (auto& point : points)
            {
                if (min_x > point.x)
                    min_x = point.x;
                if (min_y > point.y)
                    min_y = point.y;
                if (max_x < point.x)
                    max_x = point.x;
                if (max_y < point.y)
                    max_y = point.y;
            }

            diam << max_x - min_x, max_y - min_y;

            center << (min_x + max_x) / 2, (min_y + max_y) / 2, 1;
        }

        void add_point(c_point& pnt)
        {
            points.push_back(pnt);
        }

        void reset()
        {
            points.clear();
        }

        vector <c_point> points;
        double min_x, max_x;
        double min_y, max_y;
        Vector2d diam;
        Vector3d center;   // homogeneous
    };

    class c_sensor_region
    {
    public:
        double id;
        vector <double> x;
    };

    class c_world_region
    {
    public:
        double id;
        vector <double> X;
    };

    class c_actuator_command
    {
    public:
        vector<double> tau;
    };

    class c_time_slot
    {
    public:

        vector<MatrixXd> rgb_channels;
        vector <c_sensor_region> sensor_regions;
        vector <c_world_region> world_regions;
        vector <c_actuator_command> actuator_commands;

        vector <c_statistic_channel> statistic_channels;

        std::map<c_range_key, c_point_set, c_range_key> point_sets;
    };

    class c_learner
    {
    public:
        c_learner(param_context& _pc)
        {
            pc = _pc;
        }

        void read_sensors()
        {
            // output: add to the back of rgb_channels
        }

        void send_actuator_commands()
        {

        }

        void calculate_statistic_point_sets(vector <c_statistic_channel>& statistic_channels, std::map<c_range_key, c_point_set, c_range_key>& point_set_map)
        {
            int x, y;

            vector<int> ranges;

            ranges.resize(statistic_channels.size());

            point_set_map.clear();

            int width = statistic_channels[0].statistic_val.cols();
            int height = statistic_channels[0].statistic_val.rows();

            c_point pnt;
            for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                    pnt.x = x;
                    pnt.y = y;

                    for (int c = 0; c < statistic_channels.size(); c++)
                    {
                        if (statistic_channels[c].range_length == 0)
                            ranges[c] = 0;
                        else
                        {
                            ranges[c] = (statistic_channels[c].statistic_val(y, x) - statistic_channels[c].minCoeff) / statistic_channels[c].range_length;
                            if (ranges[c] < 0)  ranges[c] = 0;
                            if (ranges[c] > pc.num_chan_ranges)  ranges[c] = pc.num_chan_ranges;
                        }
                    }

                    c_range_key range_key;
                    range_key.ranges = ranges;

                    // add point to the neighbouring ranges too
                    for (int c = 0; c < statistic_channels.size(); c++)
                    {
                        // check neighbouring ranges for cases when statistic is on a range boundary
                        for (int j = -1; j <= 1; j++)
                        {
                            if (c > 0 && j == 0)
                                continue;  // otherwise not perturbed ranges is added multiple times to the same point set (num_statistic_channels to be exact)
                            c_range_key rk;
                            rk.ranges = ranges;
                            rk.ranges[c] += j;

                            c_point_set& ps = point_set_map[rk];
                            ps.add_point(pnt);
                        }
                    }
                }
            }

            for (auto it = point_set_map.begin(); it != point_set_map.end(); ++it)
            {
                it->second.calculate_boundaries();
            }
        }

        void calculate_sensor_identities()
        {
            // input: rgb_channels
            // output: sensor_regions
            auto& statistic_channels = time_slots.back().statistic_channels;
            for (int c = 0; c < time_slots.back().rgb_channels.size(); c++) {
                for (int s = 0; s < pc.num_stat_components; s++)
                {
                    statistic_channels.push_back(c_statistic_channel());
                    statistic_channels.back().calculate_image_statistic(pc, time_slots.back().rgb_channels[c], statistics[s]);
                }
            }

            calculate_statistic_point_sets(statistic_channels, time_slots.back().point_sets);
        }

        void predict_transformation_of_camera_ref_frame_from_actuator_commands()
        {
            // output: R,T with respect to the initial camera reference frame
        }

        void calculate_non_observed_world_from_behaviour()
        {

        }

        void reconcile_observed_non_observed_worlds(bool& error_significant)
        {

        }

        void calculate_observed_world()
        {
            // input: sensor_regions
            // output world_regions
            predict_transformation_of_camera_ref_frame_from_actuator_commands();
            // X ~ R1*K_inv*x1+T1; X ~ R2*(K_inv*x2)+T2 for x1,x2 of the same identity
        }

        void try_decrease_error()
        {

        }

        /*
        time
        t    sensor_input
             actuator_command
        t+1  sensor_input
             actuator_command
             reconstruction: sensor_input(cur_t-1), sensor_input, actuator_command(cur_t-1) -> observed_world_id_X
             observed_world_id_X -> non_observed_world_id_X
        t+2  sensor_input
             actuator_command
             reconstruction: sensor_input(cur_t-1), sensor_input, actuator_command(cur_t-1) -> observed_world_id_X
             behaviour: non_observed_world_id_X(cur_t-1) -> non_observed_world_id_X
             compare_world: observed_world_id_X, non_observed_world_id_X -> error, non_observed_world_id_X
             if error is significant, then modify reconstruction and/or behaviour
        t+3  sensor_input
             actuator_command
             reconstruction: sensor_input(cur_t-1), sensor_input, actuator_command(cur_t-1) -> observed_world_id_X
             behaviour: non_observed_world_id_X(cur_t-1) -> non_observed_world_id_X
             compare_world: observed_world_id_X, non_observed_world_id_X -> error, non_observed_world_id_X
             if error is significant, then modify reconstruction and/or behaviour
             ---------------------------------------------------------------------------------------------------
             The permanent assignment of X to id for observed and non_observed worlds would be a zero error solution,
             if only compare_world is used, even though a live system would perish.
             A real test should be comparing the predicted and actual sensor inputs.
        */
        void learn()
        {
            bool error_significant;

            calculate_sensor_identities();
            calculate_observed_world();
            calculate_non_observed_world_from_behaviour();
            reconcile_observed_non_observed_worlds(error_significant);
            if (error_significant)
            {
                try_decrease_error();
            }
        }

        void initialize_learning_procedures_and_parameters()
        {
            camera_params = MatrixXd::Identity(3, 3);

            statistics.clear();
            for (int i = 0; i < pc.num_stat_components; i++)
            {
                c_statistic statistic;
                statistic.calculate_linear_statistic_parameters(pc);
                statistics.push_back(statistic);
            }
        }

        void run()
        {
            initialize_learning_procedures_and_parameters();
            while (true)
            {
                time_slots.push_back(c_time_slot());
                read_sensors();
                learn();
                send_actuator_commands();
            }
        }

        param_context pc;
        list<c_time_slot> time_slots;

        MatrixXd camera_params;
        vector <c_statistic> statistics;
    };
}