/*Script to simulate random walk
 */
#include "param.h"
#include "walk.h"

int run_uniform() {
	LatticeUniform latt = LatticeUniform(3);
	std::cout << (Lattice)latt;

	WalkFixTime walk = WalkFixTime(latt);
	walk.RunWalk(
		1000, //Sample times
		100, //Time interval to sample
		1000 //Number of trajectories
	);
	walk.save("test1000-square.dat", [](double x) {return x*x;});
	return 0;
}

int run_gttest() {
	std::fstream fs("3dtest.in", std::fstream::in);
	if (!fs.good()) {
		throw std::runtime_error("Cannot open input file");
	}
	LatticeGT latt = LatticeGT();
	fs >> latt;
	fs.close();

	WalkFixTime walk = WalkFixTime(latt);

	latt.SetTemperature(300);
	std::cout << (Lattice)latt;

	walk.RunWalk(
		1000, //Sample times
		1000, //Time interval to sample
		1000 //Number of trajectories
	);
	walk.save("test1000-3dtest.dat", [](double x) {return x*x;});
	return 0;
}

int run_single(){
	std::fstream fs("bivo4-16.in", std::fstream::in);
	if (!fs.good()) {
		throw std::runtime_error("Cannot open input file");
	}
	LatticeGT latt = LatticeGT();
	fs >> latt;
	fs.close();
	
    WalkVariableTime walk = WalkVariableTime(latt);

	latt.SetTemperature(300);
	std::cout << (Lattice)latt;

    walk.RunWalk(
		100000, //Time interval to sample
		100, //Sample times
		10000 //Number of trajectories
	);
	walk.save("bivo4-16.dat", [](double x){return x*x;} );
    return 0;
}

int main(int argc, char* argv[]){
    if (argc != 2 && argc != 3){
        std::cout << "Usage : RandomWalk [-v] $PREFIX" << std::endl;
        std::cout << "Read $PREFIX.latt and $PREFIX.in as the latt file / input file and write a series of $PREFIX-*.dat" << std::endl;
        std::cout << "Use -v to display progress (must be first argument)" << std::endl;
        std::cout << "Use -l to display path selection (must be first argument)" << std::endl;
        return 0;
    }

    bool is_log = false;
    
    if (argc == 3){
        bool verbose = argv[1][0] == '-' && argv[1][1] == 'v';
        ProgressBar::is_disable = !verbose;
        is_log = argv[1][0] == '-' && argv[1][1] == 'l';
    }

    std::string prefix(argv[argc-1]);
	std::fstream fs(prefix + ".latt", std::fstream::in);
	if (!fs.good()) {
		throw std::runtime_error("Cannot open input file : latt");
	}
	LatticeGT latt = LatticeGT();
	fs >> latt;
	fs.close();

    //Check lattice
//  for (const auto& grid : latt.grids){
//      std::cout << grid.name << std::endl;
//      for (const auto& path : grid.paths){
//          std::cout <<  path.index_grid  << " " << path.shiftLattice.size() << std::endl;
//      }
//  }

    InputParam param;
	fs.open(prefix + ".in", std::fstream::in);
	if (!fs.good()) {
		throw std::runtime_error("Cannot open input file : in");
	}
    fs >> param;
    fs.close();
	
    WalkVariableTime walk = WalkVariableTime(latt);
    walk.is_log_path_select = is_log;

    walk.RunWalkAtTemperatureAdaptive(
        prefix,
        param.temperatures,
		param.num_step_per_sample, //Steps per sample
		param.num_sample, //Sample times
	    param.num_trajectory //Number of trajectories
	);
    return 0;
}
