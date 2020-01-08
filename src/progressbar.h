#pragma once
#include <iostream>
#include <string>


class ProgressBar {
    public:
        //Width : the bar length (not including name and numbers
        ProgressBar(std::string const& name, size_t max, size_t width);
        void Update(size_t val);
        void Add(size_t n);
        void Clear();
        static bool is_disable; //Global switch to disable all progressbar output

    private:
        std::string name;
        size_t max;
        size_t width;
        size_t last; //Current displayed value
};
