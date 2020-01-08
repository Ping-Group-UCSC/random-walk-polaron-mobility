#pragma once
/*
 * Control input parameters
 */
#include<iostream>
#include<vector>

class InputParam{
    public:
        int num_step_per_sample = 10;
        int num_sample = 100;
        int num_trajectory = 1000;
        std::vector<double> temperatures = {200, 250, 300, 350, 400};
        friend std::istream& operator >> (std::istream& is,  InputParam& param);
};


