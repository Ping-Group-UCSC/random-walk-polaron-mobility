#pragma once
/*
 * Define the path inside one lattice
 */
#include<iostream>
#include<fstream>
#include<string>
#include "math2.h"
#include "matrix.h"

class Grid;
/*
 * Connection between two grids
 * Including speed 
 */
struct GridPath {
    public:
        int index_grid; //Destination
        std::string name; //Name of the path (to distinguish different paths between same kinds of atoms)
        double k; //Speed constant, this will be evaluated for different temperatures
        std::vector<int> shiftLattice; // The shift of the lattice this grid goes
        double prefactor; // A prefactor for the final speed constant
//      GridPath(const Grid* grid, const std::vector& shiftLattice, double k);
		GridPath();
        GridPath(int grid, std::vector<int>&& shiftLattice, double k);
};

/*
 * One grid in the unitcell, include connections
 */
class Grid {
    public:
        int index; //Index of grid in unitcell
		std::string name; // Name of the grid
		std::string type;  //Grid type (generally for chemical enviroment)
        std::vector<double> position;
        std::vector<GridPath> paths;
        double total_rate; // Summation of rates of all possible paths
        Grid();
        Grid(int index);
        Grid(int index, const std::vector<double>& position);
};

/*
 * Represent a grid in periodical lattice
 */
struct LatticeGrid{
    public:
        std::vector<int> ix_lattice;
        int ix_grid;
        LatticeGrid(const std::vector<int>& ix_lattice, int ix_grid); 
        LatticeGrid(std::vector<int>&& ix_lattice, int ix_grid); 
};

/*
 * Represent a lattice (just graph ; no uniform assumption)
 */
class Lattice {
public:
	std::vector<Grid> grids;
	int dimension() const { return _dimension; }
	virtual std::vector<double> CalcDistance(const LatticeGrid& pos1, const LatticeGrid& pos2) const { return std::vector<double>(_dimension); } //Calculate the distance between point2 and point1 
	virtual void SetTemperature(double t) { return ; } //Set temperature for all paths
    ~Lattice();
    protected:
        int _dimension;
};
std::ostream& operator << (std::ostream& os, const Lattice & latt);


/*
 * A lattice with uniform grid points
 */
class LatticeUniform : public Lattice {
    public:
        LatticeUniform(int dim); //Create a square grid with constant speed between all nearest neighbour sites
        std::vector<int> dimensions; //Dimension of the lattice (number of grids); only used for a uniform system
        ~LatticeUniform();
        std::vector<double> CalcDistance(const LatticeGrid& pos1, const LatticeGrid& pos2) const; //Calculate the distance between point2 and point1 
		void SetTemperature(double t); //Speed are always constant ; not related to temperature;
    private:
        std::vector<int> dimension_allsize;//
        void InitDimensionSize(); // Convert dimensions to how many elements in the dimension
        int GetPeriodicalIndexIn1D(const std::vector<int> & pos); //Get the grid at given position for a uniform grid
        std::vector<int> ShiftPeriodic(const std::vector<int>& pos, const std::vector<int> &shift); //Get the shift of lattice when a grid is shift in a uniform grid
};
