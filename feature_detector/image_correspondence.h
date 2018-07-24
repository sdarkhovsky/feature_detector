#pragma once

#include "lpngwrapper.hpp"
#include "param_context.h"

#include <vector>
#include <map>
#include <random>

using namespace std;

class c_point
{
public:
    int x, y;
    double statistic[3];
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

        diam_x = max_x - min_x;
        diam_y = max_y - min_y;

        center_x = (min_x + max_x) / 2;
        center_y = (min_y + max_y) / 2;
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
    double diam_x, diam_y;
    double center_x, center_y;
};

class c_point_set_correspondence
{
public:
    int statistic_range[3];
    c_point_set* point_set1;
    c_point_set* point_set2;
};

#if 0
typedef pair<int, int> addr_type;
typedef double value_type;
typedef vector <value_type> point_value_type;


typedef void(*statistic_proc_type)(vector<c_point_set*> input, addr_type addr, c_point_set& output/*, parameters, proc_type*/);

void linear_statistic(vector<c_point_set*> input, addr_type addr, c_point_set& output)
{
    if (input.size() == 0)
        return;
    c_point_set* ps = input[0];

    const std::map<addr_type, value_type> mask = {}
}

class c_node
{
public:

    void reset()
    {
        point_set.reset();
        in_nodes.clear();
        out_nodes.clear();
    }

    c_point at_point; // at which point of an in_nodes the statistic is calculated
    statistic_proc_type statistic_proc;

    vector <double> constraint_min;
    vector <double> constraint_max;
    vector<c_node*> in_nodes;
    vector<c_node*> out_nodes;
    c_point_set point_set;
};


class c_image_graph
{
public:
    ~c_image_graph()
    {
        reset();
    }

    void reset()
    {
        for (auto& p : nodes)
        {
            if (p)
            {
                delete p;
                p = 0;
            }
        }
        nodes.clear();
    }

    void init(MatrixXd rgb_channels[], int num_channels)
    {
        reset();
        nodes.push_back(new c_node());
        nodes[0]->point_set.add_points(rgb_channels, num_channels);
    }

    void expand_graph(vector <c_node*> expanded_nodes)
    {
        expanded_nodes.clear();
        // select a statistic, statistic constraint ranges and input nodes

        // for now select the first leaf node and use a linear function sum(w[i]*v[i]) at a point's neighbourhood
        statistic_proc_type statistic_proc = linear_statistic;
        c_node* sel_node = 0;
        for (auto& p : nodes)
        {
            if (p && p->out_nodes.size() == 0)
            {
                sel_node = p;
                break;
            }
        }

        if (sel_node = 0)
            return;

        // calculate the statistic at each point of the select node's set
        for (auto& pnt : sel_node->point_set.points)
        {

        }





            // calculate the statistic, apply constraint, create new nodes
            11111111111111111111111
                c_point at_point; // at which point of an in_nodes the statistic is calculated
            void(*statistic)(vector<c_point_set> input, c_point_set& output/*, parameters, proc_type*/);
            vector <double> constraint_min;
        vector <double> constraint_max;
        111111111111111111111111
    }


    // nodes[0] points to the root node
    vector <c_node*> nodes;
};
#endif
class c_image_correspondence
{
public:
    c_image_correspondence(param_context& _pc)
    {
        pc = _pc;

        read_png_file(pc.image_path[0].c_str(), rgb_channels_1, pc.num_channels);
        read_png_file(pc.image_path[1].c_str(), rgb_channels_2, pc.num_channels);

        width = rgb_channels_1[0].cols();
        height = rgb_channels_1[0].rows();

#if 0
        image1_graph.init(rgb_channels_1, pc.num_channels);
        image2_graph.init(rgb_channels_2, pc.num_channels);
#endif
    }

#if 0
    void replicate_graph_expansion(vector <c_node*> in_expanded_nodes, c_image_graph& graph, vector <c_node*> out_expanded_nodes)
    {

    }
#endif

    void calculate_linear_statistic_parameters(int wx, int wy, MatrixXd& params)
    {
        std::default_random_engine generator;
        std::uniform_real_distribution<double> distribution(-1.0, 1.0);
        int wsize = (2 * wx + 1)*(2 * wy + 1);

        params = MatrixXd::Zero(2*wy+1, 2*wx+1);
        double sum = 0.0;
        for (int i1 = -wy; i1 <= wy; i1++)
        {
            for (int i2 = -wx; i2 <= wx; i2++)
            {
                if (i1*i1 + i2*i2 <= wsize*wsize)
                {
                    params(i1, i2) = distribution(generator);
                    if (i1 != 0 || i2 != 0)
                        sum += params(i1, i2);
                }
            }
        }

        // make total sum 0
        params(0, 0) = - sum;
    }

    void calculate_image_statistic(const MatrixXd& rgb_channel, const MatrixXd& linear_statistic_parameters, MatrixXd& statistic)
    {
        int x, y, x0, y0;
        double sum;

        int wx = (linear_statistic_parameters.cols() - 1) / 2;
        int wy = (linear_statistic_parameters.rows() - 1) / 2;

        statistic = MatrixXd::Zero(height, width);

        for (y0 = wy; y0 < height - wy; y0++) {
            for (x0 = wx; x0 < width - wx; x0++) {
                sum = 0.0;
                for (y = y0 - wy; y <= y0 + wy; y++) {
                    for (x = x0 - wx; x <= x0 + wx; x++) {
                        sum += rgb_channel(y, x)*linear_statistic_parameters(y - y0 + wy, x - x0 + wx);
                    }
                }
                statistic(y0, x0) = sum;
            }
        }
    }

    void calculate_statistic_point_sets(const MatrixXd statistic_channels[], std::map<int[3], c_point_set>& point_set_map)
    {
        int x, y;

        int range[3];
        c_point pnt;
        for (y = 0; y < height; y++) {
            for (x = 0; x < width; x++) {
                for (int c = 0; c < pc.num_channels; c++) {
                    pnt.statistic[c] = statistic_channels[c](y, x);
                    pnt.x = x;
                    pnt.y = y;
                    if (range_length[c] == 0)
                        range[c] = 0;
                    else
                    {
                        range[c] = (statistic_channels[c](y, x) - minCoeff_1[c]) / range_length[c];
                        if (range[c] < 0)  range[c] = 0;
                        if (range[c] > pc.num_chan_ranges)  range[c] = pc.num_chan_ranges;
                    }
                }

                c_point_set& ps = point_set_map[range];
                ps.add_point(pnt);
            }
        }

        for (auto it = point_set_map.begin(); it != point_set_map.end(); ++it)
        {
            it->second.calculate_boundaries();
        }
    }

    void find_RT_transformation(std::vector <c_point_set_correspondence>& correspondences)
    {

    }

    void calculate_image_correspondence()
    {
        MatrixXd linear_statistic_parameters;
        calculate_linear_statistic_parameters(pc.wx, pc.wy, linear_statistic_parameters);

        MatrixXd statistic_channels_1[3];
        MatrixXd statistic_channels_2[3];

        for (int c = 0; c < pc.num_channels; c++) {
            calculate_image_statistic(rgb_channels_1[c], linear_statistic_parameters, statistic_channels_1[c]);
            calculate_image_statistic(rgb_channels_2[c], linear_statistic_parameters, statistic_channels_2[c]);

            minCoeff_1[c] = statistic_channels_1[c].minCoeff();
            maxCoeff_1[c] = statistic_channels_1[c].maxCoeff();
            range_length[c] = (maxCoeff_1[c] - minCoeff_1[c]) / pc.num_chan_ranges;
        }

        std::map<int[3], c_point_set> point_sets_1;
        std::map<int[3], c_point_set> point_sets_2;

        calculate_statistic_point_sets(statistic_channels_1, point_sets_1);
        calculate_statistic_point_sets(statistic_channels_2, point_sets_2);

        // filter corresponence sets by set size
        std::vector <c_point_set_correspondence> correspondences;
        int statistic_localization_x = width*pc.statistic_localization;
        int statistic_localization_y = height*pc.statistic_localization;
        for (auto it = point_sets_1.begin(); it != point_sets_1.end(); ++it)
        {
            c_point_set* ps1 = &it->second;
            c_point_set* ps2 = &point_sets_2[it->first];
            if (ps1->points.size() == 0 || ps2->points.size() == 0)
                continue;
            if (ps1->diam_x < statistic_localization_x && ps1->diam_y < statistic_localization_y &&
                ps2->diam_x < statistic_localization_x && ps2->diam_y < statistic_localization_y)
            {
                c_point_set_correspondence psc;
                psc.point_set1 = ps1;
                psc.point_set2 = ps2;
                memcpy(psc.statistic_range, it->first, sizeof(psc.statistic_range));
                correspondences.push_back(psc);
            }
        }

        // filter correspondences with RANSAC
        find_RT_transformation(correspondences);
#if 0
        while (true)
        {
            vector <c_node*> image1_expanded_nodes;
            vector <c_node*> image2_expanded_nodes;
            image1_graph.expand_graph(image1_expanded_nodes);

            if (image1_expanded_nodes.size() == 0) break;

            replicate_graph_expansion(image1_expanded_nodes, image2_graph, image2_expanded_nodes);
            if (check_correspondences(image1_expanded_nodes, image2_expanded_nodes) == correspondence_failed)
            {

            }
        }
#endif
    }

    int width;
    int height;

    param_context pc;
    MatrixXd rgb_channels_1[3];
    MatrixXd rgb_channels_2[3];

    double minCoeff_1[3], maxCoeff_1[3], range_length[3];
#if 0
    c_image_graph image1_graph;
    c_image_graph image2_graph;
#endif
 };

