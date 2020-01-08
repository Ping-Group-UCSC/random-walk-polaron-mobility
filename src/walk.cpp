#include "walk.h"


//LatticePositionTime::LatticePositionTime(const std::vector<int>& ix_lattice, int ix_grid, double time):
// ix_lattice(ix_lattice), ix_grid(ix_grid), time(time){
//}

//LatticePositionTime::LatticePositionTime(std::vector<int> && ix_lattice, int ix_grid, double time):
// position(ix_lattice, ix_grid), time(time){
//}
//
LatticePositionTime::LatticePositionTime(LatticeGrid&& position, double time):
 position(position), time(time){
}

LatticePositionTime::LatticePositionTime(const LatticeGrid& position, double time):
 position(position), time(time){
}


WalkTrajectory::WalkTrajectory(const Lattice* const latt) : latt(latt){
}

void WalkTrajectory::save(std::string filename){
    if (steps.size() == 0) {
        return;
    }
    std::fstream fs(filename, std::fstream::out);
    for (const auto& step : steps){
        std::vector<double> v = latt->CalcDistance(steps[0].position, step.position);
        fs << step.time;
        std::for_each(v.begin(), v.end(), [&](double x){fs << " " << x;});
        fs << std::endl;
    }
    fs.close();
}


Walk::Walk(Lattice& latt) : latt(&latt) {}


int Walk::RandomSelect(int n){
    return (int)floor(n * ((double)random())/(random.max()+1));
}

const GridPath& Walk::RandomSelectPath(const Grid& grid){
    double pick = grid.total_rate * random() / (random.max() + 1);
    double pick0 = pick;
    for (const auto& path : grid.paths){
        pick -= path.k;
        if (pick <= 0){
//          std::cout << path.index_grid << " " << path.name << " " << tot << " " << pick0 <<  std::endl;
            return path;
        }
    }
    //Return the last
    return grid.paths.back();
}


void Walk::RunWalkAtTemperature(const std::string& prefix, const std::vector<double>& temperatures, double temperature_ref, double time_limit_ref, double num_sample, int num_trajectory) {
    ProgressBar bar = ProgressBar("Temperature ", temperatures.size(), 40);
    bar.Update(0);
    for (auto& temperature : temperatures){
        double time_limit = std::exp(std::log(time_limit_ref) * temperature_ref / temperature);
	    latt->SetTemperature(temperature);
        RunWalk(time_limit, num_sample, num_trajectory);
        save(prefix + "-t-" + std::to_string(temperature) + ".dat", [](double x){return x*x;});
        bar.Add(1);
    }
    bar.Clear();
}


void Walk::RunWalkAtTemperatureAdaptive(const std::string& prefix, const std::vector<double>& temperatures, double num_step_per_sample, double num_sample, int num_trajectory) {
    ProgressBar bar = ProgressBar("Temperature ", temperatures.size(), 40);
    bar.Update(0);
    for (auto& temperature : temperatures){
	    latt->SetTemperature(temperature);
        if (is_log_path_select){
            ClearPathSelectCounter();
        }
        //Search for a reasonable time per sample at given temperature
        double time_limit = EstimateSampleTime(num_step_per_sample);
        RunWalk(time_limit, num_sample, num_trajectory);
        save(prefix + "-t-" + std::to_string(temperature) + ".dat", [](double x){return x*x;});
        bar.Add(1);

        if (is_log_path_select){
            PrintPathSelectCounter();
        }
    }
    bar.Clear();
}


void Walk::RunWalk(double time_limit, double num_sample, int num_trajectory) {
    this->time_limit = time_limit;
	this->num_sample = num_sample;
    this->num_trajectory = num_trajectory; 
    if (is_keep_trajectory()){
        trajectories.clear();
        trajectories.reserve(num_trajectory);
    }
    else{
        distance_temp.clear();
        distance_temp.reserve(num_trajectory);
    }

    ProgressBar bar = ProgressBar("Trajectory  ", num_trajectory, 40);
    bar.Update(0);

    for (int i = 0 ; i < num_trajectory ; i ++){
        //Save the distance
        if (!is_keep_trajectory()){
            distance_temp.emplace_back(InterpolateDistance(RunWalkOnce(time_limit, num_sample)));
        }
        else{
            trajectories.emplace_back(RunWalkOnce(time_limit, num_sample));
        }
        bar.Add(1);
    }

    bar.Clear();
}

void Walk::AddPathSelectCounter(const Grid& grid, const GridPath& path) {
    std::string name = grid.type + latt->grids[path.index_grid].type + path.name;
	if (dic_path_select_counter.find(name) == dic_path_select_counter.end()){
        dic_path_select_counter[name] = 0;
    }
    dic_path_select_counter[name] += 1;
    return;
}

void Walk::PrintPathSelectCounter() {
    for (const auto& kvp : dic_path_select_counter){
        std::cout << kvp.first << " " << kvp.second << std::endl;
    }
}

void Walk::ClearPathSelectCounter() {
   dic_path_select_counter = std::map<std::string, int>();
}

WalkFixTime::WalkFixTime(Lattice& latt) : Walk(latt) {}

void WalkFixTime::save(std::string filename, std::function<double(double)> const& func_x) {
	if (num_trajectory == 0) return;
	//Check if all trajectories have same length (not used; trajectories may have variable lengths)
	
	for (const auto& traj : trajectories) {
		if (traj.steps.size() != trajectories[0].steps.size()) {
			throw std::runtime_error("Trajectories must have same length");
		}
	}
	std::fstream fs(filename, std::fstream::out);
	int n = trajectories[0].steps.size();
	for (int i = 0; i < n; i++) {
//      std::valarray<double> avg(num_trajectory);
//  	std::transform(trajectories.begin(), trajectories.end(), std::begin(avg),
//  		[func_x, i, this](const WalkTrajectory& traj){
//  		return func_x(latt->CalcDistance(traj.steps[0].position, traj.steps[i].position));
//  	});
//  	fs << trajectories[0].steps[i].time << " " 
//          << avg.sum() / avg.size() << " " 
//          << avg.min() << " " 
//          << avg.max() << std::endl;
	}
	fs.close();
}

WalkTrajectory WalkFixTime::RunWalkOnce(double time_limit, double time_sample){
//Random seed current time
    random.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());

    std::vector<int> pos0(latt->dimension()); //Initial at lattice (0,0,0,0...)
    std::fill(pos0.begin(), pos0.end(), 0);

    //Random starting grid
    LatticePositionTime posNow(
            LatticeGrid(pos0, RandomSelect(latt->grids.size())), 
            0);
    auto traj = WalkTrajectory(latt);

    double time_now = 0;
    double time_last_sample = 0;
    const Grid* grid = &(latt->grids[posNow.position.ix_grid]);
    while (time_now < time_limit){
        //Record
        for (; time_last_sample <= time_now ; time_last_sample += time_sample ){
            traj.steps.push_back(posNow);
        }

        //Select path
		const auto& path = RandomSelectPath(*grid);
        time_now += 1/path.k;
        std::transform(
                posNow.position.ix_lattice.begin(), 
                posNow.position.ix_lattice.end(), 
                path.shiftLattice.begin(), 
                posNow.position.ix_lattice.begin(), 
                std::plus<int>());
        grid = &(latt->grids[path.index_grid]);
        posNow.position.ix_grid = grid->index;
        posNow.time = time_now;
    }
    return traj;
}

WalkVariableTime::WalkVariableTime(Lattice& latt) : Walk(latt) {}


std::vector<std::vector<double>> WalkVariableTime::InterpolateDistance(const WalkTrajectory& traj){
	double time_max = time_limit * num_sample;
	double time_step = time_limit;
    const auto& steps = traj.steps;
    //Convert to distance
    std::vector<std::vector<double>> distance_alldim(steps.size(), std::vector<double>(latt->dimension()));
    std::transform(steps.begin(), steps.end(), distance_alldim.begin(), [this, steps](auto& x) {return latt->CalcDistance(steps[0].position, x.position);});

    std::vector<std::vector<double>> avg(latt->dimension());

    //Linear interpolation
    for (int idim = 0 ; idim < latt->dimension() ; idim ++){
        avg[idim].resize(num_sample);
        for (int i = 0; i < steps.size() - 1; i++) {
            //Find all points in one time step
            for (int time = std::ceil(steps[i].time / time_step); time <= steps[i + 1].time / time_step && time < num_sample; time++) {
                double t2 = time * time_step;
                avg[idim][time] = 
                    (distance_alldim[i][idim] * (steps[i + 1].time - t2) + distance_alldim[i + 1][idim] * (t2 - steps[i].time))
                    / (steps[i+1].time - steps[i].time);
            }
        }
    }
    return avg;
}

void WalkVariableTime::save(std::string filename, std::function<double(double)> const & func_x)
{
	if (num_trajectory == 0) return;
	//Try to find a reasonable linear grid to interpolate
	double time_max = time_limit * num_sample;
	std::vector<std::vector<std::valarray<double>>> avg(num_sample);
    for (int i = 0 ; i < num_sample; i ++){
        avg[i] = std::vector<std::valarray<double>>(latt->dimension(), std::valarray<double>(num_trajectory));
    }
	double time_step = time_max / avg.size();

    double num_step_per_sample = 0;
	//For all trajectories
    for (int it = 0 ; it < num_trajectory ; it ++){
        std::vector<std::vector<double>> avgt;
        if (is_keep_trajectory()){
            num_step_per_sample += trajectories[it].steps.size();
            avgt = InterpolateDistance(trajectories[it]);
        }
        else{
            avgt = distance_temp[it];
        }
        for (int time = 0 ; time < num_sample ; time++){
            for (int idim = 0 ; idim < latt->dimension() ; idim ++){
                avg[time][idim][it] = func_x(avgt[idim][time]);
            }
        }
	}
    num_step_per_sample /= (num_trajectory * num_sample);
    if (num_step_per_sample > 0) {
        std::cout << "Steps per sample: " << num_step_per_sample << std::endl;
    }
	//Output 
	std::fstream fs(filename, std::fstream::out);

	for (int i = 0; i < avg.size(); i++) {
		fs << i * time_step << " " ;
        for (int idim = 0 ; idim < latt->dimension() ; idim ++){
            fs << avg[i][idim].sum() / num_trajectory << " " ; 
        }
        fs << std::endl;
	}
	fs.close();
}

WalkTrajectory WalkVariableTime::RunWalkOnce(double time_limit, double num_sample)
{
	//Random seed current time
	random.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());

	std::vector<int> pos0(latt->dimension()); //Initial at lattice (0,0,0,0...)
	std::fill(pos0.begin(), pos0.end(), 0);

	//Random starting grid
	LatticePositionTime posNow(
		LatticeGrid(pos0, RandomSelect(latt->grids.size())),
		0);
	auto traj = WalkTrajectory(latt);

	double time_now = 0;
	double time_last_sample = 0;
	const Grid* grid = &(latt->grids[posNow.position.ix_grid]);
	while (time_now < time_limit * num_sample ) { //Run until enough time
		//Select path
		const auto& path = RandomSelectPath(*grid);
        if (is_log_path_select){
            AddPathSelectCounter(*grid, path);
        }
		time_now += - 1 / grid->total_rate * log(((double)random()+1)/(random.max()+1));

		std::transform(
			posNow.position.ix_lattice.begin(),
			posNow.position.ix_lattice.end(),
			path.shiftLattice.begin(),
			posNow.position.ix_lattice.begin(),
			std::plus<int>());
		grid = &(latt->grids[path.index_grid]);
		posNow.position.ix_grid = grid->index;
		posNow.time = time_now;
		//Record every step near give time
		traj.steps.push_back(posNow);
	}
	return traj;
}

double WalkVariableTime::EstimateSampleTime(double num_step_per_sample)
{
	//Random seed current time
	random.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());

	std::vector<int> pos0(latt->dimension()); //Initial at lattice (0,0,0,0...)
	std::fill(pos0.begin(), pos0.end(), 0);

	//Random starting grid
	LatticePositionTime posNow(
		LatticeGrid(pos0, RandomSelect(latt->grids.size())),
		0);

//Run n_test to esimate 
    int n_test = 100;

	double time_now = 0;
	double time_last_sample = 0;
	const Grid* grid = &(latt->grids[posNow.position.ix_grid]);
    for (int i = 0 ; i < num_step_per_sample * n_test; i ++){
		//Select path
		const auto& path = RandomSelectPath(*grid);
		time_now += 1 / grid->total_rate;
		std::transform(
			posNow.position.ix_lattice.begin(),
			posNow.position.ix_lattice.end(),
			path.shiftLattice.begin(),
			posNow.position.ix_lattice.begin(),
			std::plus<int>());
		grid = &(latt->grids[path.index_grid]);
		posNow.position.ix_grid = grid->index;
		posNow.time = time_now;
		//Record every step near give time
    }
    return posNow.time / n_test; 
}



