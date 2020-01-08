#pragma once
/*
* This is for lattice with inhomogeneous grid, and paths are determined by grid type
*/

#include "lattice.h"

class PathBarrier {
public:
    PathBarrier();
    PathBarrier(double Ea, double prefactor);
    ~PathBarrier();
    double Ea; //Activation energy
    double prefactor; // Prefactor
};

/*
* A lattice with inhomogeneous grid points
* There is a fixed activation energy between two kinds of points
*/
class LatticeGT : public Lattice {
public:
	LatticeGT(); //Create a square grid with constant speed between all nearest neighbour sites
    std::vector<double> CalcDistance(const LatticeGrid& pos1, const LatticeGrid& pos2) const; //Calculate the distance between point2 and point1 
	~LatticeGT();
	void SetTemperature(double t); //Set all speed based on temperature
	friend std::istream& operator >> (std::istream& is, LatticeGT & latt);
	friend std::ostream& operator << (std::ostream& os, const LatticeGT & latt);
private:
	Matrix<double> cell;
	std::map<std::string, int> map_name_index; //AtomName - index
    std::map<std::string, std::map<std::string, std::map<std::string, PathBarrier>>> map_barrier; //AtomType starting - AtomType ending to speed
};
