#include "latticegt.h"

PathBarrier::PathBarrier() {}

PathBarrier::~PathBarrier() {}

PathBarrier::PathBarrier(double Ea, double prefactor) :
   Ea(Ea), prefactor(prefactor)
{
}

LatticeGT::LatticeGT() : cell(3,3)
{
	_dimension = 3;
}

LatticeGT::~LatticeGT()
{
}

std::vector<double> LatticeGT::CalcDistance(const LatticeGrid & pos1, const LatticeGrid & pos2) const
{
	Matrix<double> delta;
	delta = Math2::to_double(pos2.ix_lattice - pos1.ix_lattice);
	Matrix<double> d2 =  Matrix<double>::dot(cell, delta);
    std::vector<double> v;
    v.resize(dimension());
    for (int i = 0 ; i < dimension() ; i++){
        v[i] = grids[pos2.ix_grid].position[i] - grids[pos1.ix_grid].position[i] + d2(i,0);
    }
    return v;
}

std::istream & operator>>(std::istream & is, LatticeGT & latt)
{
	//Check first line
	std::string st;
	std::getline(is, st);
	if (st != "3D") {
		throw std::runtime_error("Input file does not start with '3D'");
	}
	//Lattice vector : row-wise
	for (size_t i = 0 ; i < 3 ; i ++)
		for (size_t j = 0; j < 3; j++) {
			is >> (latt.cell(j, i));
		}
	int num_site = 0;
	//Number of atoms;
	//Format : atom name ; atom type ; atom coordinates in Cartesian (3 numbers) ; number of paths ; path to name + lattice shift (3 numbers);
	//Name and type are both string
	//Name cannot be duplicated
	is >> num_site;
	if (num_site == 0) {
		throw std::runtime_error("Number of sites should not be zero");
	}
	double unit_pos;
	is >> unit_pos;
	latt.map_name_index.clear();
	latt.map_barrier.clear();
	latt.grids.resize(num_site);
	std::vector<std::vector<std::string>> pathnames(num_site);//Store the name; set the value after read all sites
	for (int i = 0; i < num_site; i++) {
		Grid& grid = latt.grids[i];
		grid.position.resize(3);
		is >> grid.name;
		is >> grid.type >> grid.position[0] >> grid.position[1] >> grid.position[2];
		std::transform(grid.position.begin(), grid.position.end(), grid.position.begin(), [unit_pos](double x) {return x * unit_pos;});
		grid.index = i;
		if (latt.map_name_index.find(grid.name) != latt.map_name_index.end()) {
			throw std::runtime_error("Duplicated name");
		}
		latt.map_name_index[grid.name] = i;
		int num_path = 0;
		is >> num_path;
		grid.paths.resize(num_path);
		for (int j = 0; j < num_path; j++) {
            auto& pathj = grid.paths[j];
			std::string st2;
			is >> st2;
			pathnames[i].push_back(st2);
			pathj.shiftLattice.resize(3);
			is >> pathj.name >>  pathj.shiftLattice[0] >> pathj.shiftLattice[1] >> pathj.shiftLattice[2];
		}
	}
	//Set paths
	for (size_t i = 0; i < num_site; i++) {
		for (size_t j = 0; j < pathnames[i].size(); j++) {
			latt.grids[i].paths[j].index_grid = latt.map_name_index[pathnames[i][j]];
		}
	}
	//Read path activation energies
	int num_path = 0;
	is >> num_path;
	if (num_path == 0) {
		throw std::runtime_error("Number of paths should not be zero");
	}
	for (int i = 0; i < num_path; i++) {
		std::string st1, st2, stp;
		double s1, s2;
        double prefactor;
		is >> st1 >> st2 >> stp >> s1 >> s2 >> prefactor;
		if (latt.map_barrier.find(st1) == latt.map_barrier.end()) {
			latt.map_barrier[st1] = std::map<std::string, std::map<std::string, PathBarrier>>();
		}
        if (latt.map_barrier[st1].find(st2) == latt.map_barrier[st1].end()){
			latt.map_barrier[st1][st2] = std::map<std::string, PathBarrier>();
        }

		if (latt.map_barrier.find(st2) == latt.map_barrier.end()) {
			latt.map_barrier[st2] = std::map<std::string, std::map<std::string, PathBarrier>>();
		}
        if (latt.map_barrier[st2].find(st1) == latt.map_barrier[st2].end()){
			latt.map_barrier[st2][st1] = std::map<std::string, PathBarrier>();
        }

		latt.map_barrier[st1][st2][stp] = PathBarrier(s1, prefactor);
		latt.map_barrier[st2][st1][stp] = PathBarrier(s2, prefactor);
	}
	return is;
}

void LatticeGT::SetTemperature(double t) {
	for (auto& grid : grids) {
        grid.total_rate = 0;
		for (auto& path : grid.paths) {
//          std::cout << grid.type << grids[path.index_grid].type << path.name << std::endl;
			PathBarrier& p = map_barrier[grid.type][grids[path.index_grid].type][path.name];
			path.k = p.prefactor * std::exp(- p.Ea / t * (1.60217662e-19 / 1.38064852e-23)); //eV
            grid.total_rate += path.k;
		}
	}
}

std::ostream & operator<<(std::ostream & os, const LatticeGT & latt)
{
	return os;
}
