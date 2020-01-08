#include "param.h"



std::istream& operator >> (std::istream& is,  InputParam& param){
    int n_temperature;
    is >> param.num_step_per_sample
        >> param.num_sample 
        >> param.num_trajectory 
        >> n_temperature;

    param.temperatures.resize(n_temperature);

    for (int i = 0 ; i < n_temperature ; i ++){
        is >> param.temperatures[i];
    }
}

