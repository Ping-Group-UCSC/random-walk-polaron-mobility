/*
 * Some simple math operations
 */
#pragma once
//#define _SCL_SECURE_NO_WARNINGS //For MSCV with std::copy; defined in project file
#pragma warning(disable:4996)

#include<vector>
#include<iterator>
#include<algorithm>
#include<cmath>
#include<queue>
#include<set>
#include<stack>
#include<list>
#include<valarray>
#include<map>

class Math2 {
    public:
        template<class T>
        static std::vector<double> to_double(const std::vector<T>& v){
            return std::vector<double>(v.begin(), v.end());
        }

        template<class T, class BinaryOperation>
        static void vector_op(const std::vector<T>& v1, 
                const std::vector<T>& v2, 
                std::vector<T>& v3, 
                BinaryOperation op){
            std::transform(v1.begin(), v1.end(),
                    v2.begin(), v3.begin(),
                    op);
        }
        template<class T, class BinaryOperation>
        static void vector_op_inplace(std::vector<T>& v1, 
                const std::vector<T>& v2, 
                BinaryOperation op){
            std::transform(v1.begin(), v1.end(),
                    v2.begin(), v1.begin(),
                    op);
        }
};

template<typename T>
std::vector<T> operator* (const std::vector<T>& v1, const std::vector<T>& v2){
    std::vector<T> v(v1.size());
    Math2::vector_op<T>(v1, v2, v, std::multiplies<T>());
    return v;
}

template<typename T>
std::vector<T> operator+ (const std::vector<T>& v1, const std::vector<T>& v2){
    std::vector<T> v(v1.size());
    Math2::vector_op<T>(v1, v2, v, std::plus<T>());
    return v;
}

template<typename T>
std::vector<T> operator- (const std::vector<T>& v1, const std::vector<T>& v2){
    std::vector<T> v(v1.size());
    Math2::vector_op<T>(v1, v2, v, std::minus<T>());
    return v;
}

template<typename T>
std::vector<T> operator/ (const std::vector<T>& v1, const std::vector<T>& v2){
    std::vector<T> v(v1.size());
    Math2::vector_op<T>(v1, v2, v, std::divides<T>());
    return v;
}
