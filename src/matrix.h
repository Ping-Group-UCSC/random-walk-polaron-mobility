#pragma once
/*
* Matrix type support simple operations
*/
#include "math2.h"

template<typename T>
class Matrix
{
public:
	Matrix() {}
	Matrix(int i, int j) :
		data(i*j), nRow(i), nCol(j) {}
	~Matrix() {};

	//Vector always initialized as a Column-vector 
	//Note there is no possible move constructor for vector -> valarray
	Matrix<T>& operator= (const std::vector<T>& v) {
		data.resize(v.size());
		std::copy(v.begin(), v.end(), std::begin(data));
		nRow = v.size();
		nCol = 1;
		return *this;
	}

	T& operator()(size_t i, size_t j) { return data[nCol*i + j]; }

	//Norm of vector (assume 1D)
	double norm() {
		return sqrt((data * data).sum());
	}

	static Matrix<T> dot(const Matrix<T>& m1, const Matrix<T>& m2) {
		if (m1.nCol != m2.nRow) {
			throw std::runtime_error("Inconsistent matrix size");
		}
		int n = m1.nCol;
		Matrix<T> m3(m1.nRow, m2.nCol);
		for (size_t i = 0 ;  i  < m1.nRow ; i ++)
			for (size_t j = 0; j < m2.nCol; j++){
				m3(i,j) = (m1.data[std::slice(i*m1.nCol, n, 1)] * m2.data[std::slice(j, n, m2.nCol)]).sum();
			}
		return m3;
	}

private:
	size_t nRow; //Dimension 1
	size_t nCol; // Dimension 2
	std::valarray<T> data;
};

