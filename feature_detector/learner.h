#pragma once

#include "lpngwrapper.hpp"
#include "param_context.h"

#include <vector>
#include <map>
#include <random>
#include <iostream>
#include <fstream>  

using namespace std;

class c_learner
{
public:
    c_learner(param_context& _pc)
    {
        pc = _pc;
    }

    void read_sensors()
    {

    }

    void send_commands()
    {

    }

    void learn()
    {

    }

    void run()
    {
        while (true)
        {
            read_sensors();
            send_commands();
            learn();
        }
    }

    param_context pc;
};
