/*
 * Define the path inside one lattice
 */
#include "lattice.h"

/*
GridPath::GridPath(const Grid * grid, const std::vector<int>& shiftLattice, double k)
    shiftLattice(shiftLattice0)
{
    this->grid = grid;
    this->k = k;
}
*/

GridPath::GridPath()
{
}

GridPath::GridPath(int grid, std::vector<int>&& shiftLattice0, double k) :
    shiftLattice(shiftLattice0), index_grid(grid), k(k), prefactor(1)
{
}

Grid::Grid() : index(0){
}

Grid::Grid(int index) : index(index){
}

Grid::Grid(int index, const std::vector<double>& position) : 
    index(index), position(position)
{}


LatticeGrid::LatticeGrid(std::vector<int>&& ix_lattice, int ix_grid):
    ix_lattice(ix_lattice), ix_grid(ix_grid) {}

LatticeGrid::LatticeGrid(const std::vector<int>& ix_lattice, int ix_grid):
    ix_lattice(ix_lattice), ix_grid(ix_grid) {}

Lattice::~Lattice() {}

LatticeUniform::~LatticeUniform() {}

std::vector<double> LatticeUniform::CalcDistance(const LatticeGrid& pos1, const LatticeGrid& pos2) const{
    std::vector<double> delta(dimension());
    delta = Math2::to_double((pos2.ix_lattice - pos1.ix_lattice)*dimensions) + (grids[pos2.ix_grid].position - grids[pos1.ix_grid].position);
    return delta;
}

void LatticeUniform::SetTemperature(double t)
{
	for (auto& grid : grids) {
		for (auto& path : grid.paths) {
			path.k = 1;
		}
	}
}

void LatticeUniform::InitDimensionSize(){
    int a = 1;
    dimension_allsize.resize(dimension());
    for (int i = dimension()-1 ; i >=0 ; i --){
        dimension_allsize[i] = a;
        a *= dimensions[i];
    }
    return;
}

std::vector<int> LatticeUniform::ShiftPeriodic(const std::vector<int>& pos, const std::vector<int> &shiftGrid){
    //Check prerequisite
    if (dimension() != pos.size() || dimension() != shiftGrid.size()){
        throw std::runtime_error("Size mismatch");
    }
    std::vector<int> shiftLattice(dimension());
    for (size_t i = 0 ; i < pos.size() ; i ++){
        int p2 = pos[i] + shiftGrid[i];
        int p3 = (p2  % dimensions[i] + dimensions[i]) % dimensions[i];
        shiftLattice[i] = (p3-p2) / dimensions[i];
        if ( (p3 - p2)  % dimensions[i] != 0){
            throw std::runtime_error("Division error");
        }
    }
    return shiftLattice;
}

int LatticeUniform::GetPeriodicalIndexIn1D(const std::vector<int> & pos){
    //Check prerequisite
    if (dimension() != pos.size()){
        throw std::runtime_error("Size mismatch");
    }
    if (dimension_allsize.size() == 0){
        InitDimensionSize();
    }
    int a = 0;
    for (size_t i = 0 ; i < pos.size() ; i++){
//Use periodical condition
//This (a%b+b)%b ensures the result in 0-b (otherwise relies on sign)
        int p = (pos[i] % dimensions[i] + dimensions[i]) % dimensions[i];
        a += p * dimension_allsize[i];
    }
    return a;
}


LatticeUniform::LatticeUniform(int dim) {
	int size = 2; //We test 2 for 1D, 2x2 for 2D and 2x2x2 for 3D ...
	_dimension = dim;
	dimensions.resize(dim);
	std::fill(dimensions.begin(), dimensions.end(), size);
	grids.resize((size_t)std::pow(size, dim));
	//Initialize cell parameters


	//Generate all possible nearest neighbours
	std::vector<std::vector<int>> neighbours(0);

	for (int i = 0; i < dim; i++) {
		std::vector<int> shift(dim);
		std::fill(shift.begin(), shift.end(), 0);
		for (auto j : std::vector<int>({ -1,1 })) {
			shift[i] = j;
			neighbours.push_back(shift);
		}
	}

	//Generate all grid points
	//Generate all indicies
	std::list<std::vector<int>> q_index;
	q_index.emplace_back(std::vector<int>({}));
	for (int n = 0; n < dim; n++) {
		int n2 = q_index.size();
		for (int j = 0; j < n2; j++) {
			std::vector<int> vec = q_index.front();
			q_index.pop_front();
			vec.push_back(0);
			for (int i = 0; i < size; i++) {
				vec[vec.size() - 1] = i;
				q_index.push_back(vec);
			}
		}
	}

	for (auto& index_all : q_index)
	{
        int ix = GetPeriodicalIndexIn1D(index_all);
        auto& grid = grids[ix];
        grid.index = ix;
        grid.position = std::vector<double>(index_all.begin(), index_all.end());
        grid.paths.clear();
        grid.paths.reserve(neighbours.size()); //2D system always have 4 neighbours
        for (auto shift : neighbours){
            grid.paths.emplace_back(
                        GridPath(
                        GetPeriodicalIndexIn1D(index_all + shift),
                        ShiftPeriodic(index_all, shift),
                        1
                        ));
        }
    }
}

std::ostream& operator << (std::ostream& os, const Lattice & latt){
    os << "Number of grids : " << latt.grids.size() << std::endl;
    for (auto grid : latt.grids){
        os << "Grid " << grid.index << std::endl;
        for (auto path : grid.paths){
            os << "    To " << latt.grids[path.index_grid].index <<  " Speed "  << path.k << std::endl;
        }
    }
    return os;
}

