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
    static const int num_arm_joints = 3;
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

    class c_actuator
    {
    public:
        c_actuator()
        {
            for (int i = 0; i < num_arm_joints; i++)
            {
                joint_RT[i] = Matrix4d::Identity();
            }

            joint_axis[0] << 0, 1, 0;
            for (int i = 1; i < num_arm_joints; i++)
            {
                joint_axis[i] << 1, 0, 0;
            }

            for (int i = 0; i < num_arm_joints; i++)
            {
                joint_orig[i] << 0, 0, 0;
            }

            ang_mult = 1.0;

            calculate_RT_from_joint_RT();
        }

        void calculate_joint_transformations(const c_actuator& prev_actuator)
        {
            auto& prev_ang_vel = prev_actuator.ang_vel;
            double dt = 1.0;
            for (int i = 0; i < num_arm_joints; i++)
            {
                Matrix4d M;
                M << 1, -dt*prev_ang_vel[i]*joint_axis[i](2), dt*prev_ang_vel[i] * joint_axis[i](1), joint_orig[i](0),
                    dt*prev_ang_vel[i] * joint_axis[i](2), 1, -dt*prev_ang_vel[i] * joint_axis[i](0), joint_orig[i](1),
                    -dt*prev_ang_vel[i] * joint_axis[i](1), dt*prev_ang_vel[i] * joint_axis[i](0), 1, joint_orig[i](2),
                    0, 0, 0, 1;

                joint_RT[i] = M*prev_actuator.joint_RT[i];  // (NB 15, p.136)
            }

            calculate_RT_from_joint_RT();
        }

        void get_RT_inv(Matrix3d& R, Vector3d& T)
        {
            R = RT_inv.topLeftCorner(3, 3);
            T = RT_inv.topRightCorner(3, 1);
        }

        Vector3d joint_axis[num_arm_joints];
        Vector3d joint_orig[num_arm_joints];
        MatrixXd joint_RT[num_arm_joints];   // transformation from a reference frame fixed to the link k to a reference frame fixed to the link k+1
        double ang_vel[num_arm_joints];
        double ang_mult;
        MatrixXd RT; // transformation from a reference frame fixed to the link 0 to a reference frame fixed to the link num_arm_joints
        MatrixXd RT_inv;

    protected:
        void calculate_RT_from_joint_RT()
        {
            RT = joint_RT[0];
            for (int i = 1; i < num_arm_joints; i++)
            {
                RT = joint_RT[i] * RT;
            }
            RT_inv = RT.inverse();
        }
    };

    class c_time_slot
    {
    public:
        c_time_slot()
        {
            K_inv = MatrixXd::Identity(3, 3);
        }

        vector<MatrixXd> rgb_channels;
        vector <c_sensor_region> sensor_regions;
        vector <c_world_region> world_regions;

        vector <c_statistic_channel> statistic_channels;

        std::map<c_range_key, c_point_set, c_range_key> point_sets;
        std::map<c_range_key, Vector3d, c_range_key> obs_X;
        std::map<c_range_key, Vector3d, c_range_key> pred_X;

        c_actuator actuator;
        MatrixXd K_inv;
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
            std::uniform_int_distribution<> distribution(-1, 1);
            if (pc.seed_random_engines)
            {
                actuator_command_generator.seed(time(NULL));
            }

            auto cur_ts = time_slots.rbegin();
            for (int i = 0; i < num_arm_joints; i++)
            {
                cur_ts->actuator.ang_vel[i] = distribution(actuator_command_generator);
            }
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

            int width = time_slots.back().rgb_channels[0].cols();
            int height = time_slots.back().rgb_channels[0].rows();

            auto& statistic_channels = time_slots.back().statistic_channels;
            for (int c = 0; c < time_slots.back().rgb_channels.size(); c++) {
                for (int s = 0; s < pc.num_stat_components; s++)
                {
                    statistic_channels.push_back(c_statistic_channel());
                    statistic_channels.back().calculate_image_statistic(pc, time_slots.back().rgb_channels[c], statistics[s]);
                }
            }

            calculate_statistic_point_sets(statistic_channels, time_slots.back().point_sets);

            // filter corresponence sets by set size
            int statistic_localization_x = width*pc.statistic_localization;
            int statistic_localization_y = height*pc.statistic_localization;
            auto& point_sets = time_slots.back().point_sets;
            for (auto it = point_sets.begin(); it != point_sets.end(); )
            {
                const c_range_key* prk = &it->first;
                c_point_set* ps = &it->second;

                if (ps->points.size() == 0)
                {
                    it = point_sets.erase(it);
                    continue;
                }

                if (!(ps->diam(0) < statistic_localization_x && ps->diam(1) < statistic_localization_y))
                {
                    it = point_sets.erase(it);
                    continue;
                }
                ++it;
            }
        }

        void calculate_non_observed_world_from_behaviour()
        {
            auto cur_ts = time_slots.rbegin();
            auto prev_ts = std::next(cur_ts);
            cur_ts->pred_X = prev_ts->obs_X;
        }

        void reconcile_observed_non_observed_worlds(double& max_diff)
        {
            auto cur_ts = time_slots.rbegin();
            max_diff = 0;
            for (auto it = cur_ts->obs_X.begin(); it != cur_ts->obs_X.end(); ++it)
            {
                const c_range_key* prk = &it->first;
                auto& obs_X = it->second;

                auto& pred_X = cur_ts->pred_X[*prk];
                double diff = (obs_X - pred_X).norm();

                if (diff > max_diff)
                    max_diff = diff;
            }
        }

        void calculate_observed_world()
        {
            // input: sensor_regions
            // output world_regions
            auto cur_ts = time_slots.rbegin();
            auto prev_ts = std::next(cur_ts);

            cur_ts->actuator.calculate_joint_transformations(prev_ts->actuator);

            for (auto it = prev_ts->point_sets.begin(); it != prev_ts->point_sets.end(); ++it)
            {
                const c_range_key* prk = &it->first;
                c_point_set* prev_ps = &it->second;

                c_point_set* cur_ps = &cur_ts->point_sets[*prk];

                if (prev_ps->points.size() == 0 || cur_ps->points.size() == 0)
                    continue;
                // calculate obs_X in the reference frame of the camera at previous time slot (NB 15,pp.133, 136)
                Matrix3d R;
                Vector3d T;

                prev_ts->actuator.get_RT_inv(R, T);
                Vector3d p1 = T;
                Vector3d d1 = R*prev_ts->K_inv*prev_ps->center;

                cur_ts->actuator.get_RT_inv(R, T);
                Vector3d p2 = T;
                Vector3d d2 = R*cur_ts->K_inv*cur_ps->center;
                Vector3d n2 = d2.cross(d1.cross(d2));
                cur_ts->obs_X[*prk] = p1 + (p2 - p1).dot(n2) / (d1.dot(n2))*d1;
            }
        }

        void try_decrease_error()
        {
            /*
            K_inv, ang_mult, joint_axis, joint_orig,
            joint_RT[i] for time slot 0
            statistics
            */
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
            double max_diff;

            calculate_sensor_identities();
            if (time_slots.size() <= 1)
                return;

            calculate_observed_world();
            calculate_non_observed_world_from_behaviour();
            reconcile_observed_non_observed_worlds(max_diff);
            try_decrease_error();
        }

        void initialize_learning_procedures_and_parameters()
        {
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

        vector <c_statistic> statistics;

        std::default_random_engine actuator_command_generator;
    };
}