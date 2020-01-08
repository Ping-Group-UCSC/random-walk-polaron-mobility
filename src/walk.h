#pragma once
/*
 * Random walk 
 */
#include "lattice.h"
#include "latticegt.h"
#include "progressbar.h"
#include<random>
#include<chrono>
#include<functional>

struct LatticePositionTime {
    public:
        LatticeGrid position;
        double time;
//      LatticePositionTime(const std::vector<int>& indexLattice, int indexGrid, double time);
        LatticePositionTime(const LatticeGrid& position, double time);
        LatticePositionTime(LatticeGrid&& position, double time);
};

class WalkTrajectory {
    public:
        const Lattice* const latt;
        std::vector<LatticePositionTime> steps;
        void save(std::string filename); //Save to given file
        WalkTrajectory(const Lattice* const latt);
		
};

class Walk {
    public:
        Walk(Lattice& latt0);
        //Run lots of simulations at different temperatures
        //A temperature_ref is provided to setup time_limit for different temperatures (as it may vary a lot due to exponential relationship
        void RunWalkAtTemperature(const std::string& prefix, const std::vector<double>& temperatures, double temperature_ref, double time_limit_ref, double num_sample, int num_trajectory);
        void RunWalkAtTemperatureAdaptive(const std::string& prefix, const std::vector<double>& temperatures, double num_step_per_sample, double num_sample, int num_trajectory);
        void RunWalk(double time_limit, double num_sample, int num_trajectory);
        virtual WalkTrajectory RunWalkOnce(double time_limit, double num_sample) {throw std::runtime_error("Call dummy function");}
		virtual void save(std::string filename, std::function<double(double)> const& func_x) {throw std::runtime_error("Call dummy function");}
        int RandomSelect(int n); //Randomly select one from given number
        const GridPath& RandomSelectPath(const Grid& grid); //Randomly select one path based on 1/k
        std::mt19937 random;
        bool is_log_path_select; //If we log the path selected
    protected:
        Lattice* const latt;
		double num_sample;
        double time_limit;
        int num_trajectory;
        std::vector<WalkTrajectory> trajectories; // Record all trajectories


        std::vector<std::vector<std::vector<double>>> distance_temp; //Store the interpolated trajectory distance
		virtual std::vector<std::vector<double>> InterpolateDistance(const WalkTrajectory& traj) { throw std::runtime_error("Not implemented"); } //Interpolate trajectory to give time grid
        virtual bool is_keep_trajectory() {return true;} // Stroe /not store the full trajectory; just store the sample points
        virtual double EstimateSampleTime(double num_step_per_sample) {throw std::runtime_error("Not implemented");} //Estimate the temperature used to estimate the sample

        void PrintPathSelectCounter(); //Print number of times to select a specific path
        void ClearPathSelectCounter(); //Clear the counter
        void AddPathSelectCounter(const Grid& grid, const GridPath& path); //Log the path selected
        std::map<std::string, int> dic_path_select_counter = std::map<std::string, int>(); //Record number of times to select a specific path
};

//Deal with fixed timestep (all speed constant are the same in lattice)
class WalkFixTime : public Walk {
    public:
        WalkFixTime(Lattice& latt0);
		//Save average of all trajectories assuming time grid in all trajectories are the same; a function convert distance to something is passed and averaged
        WalkTrajectory RunWalkOnce(double time_limit, double num_sample);
		void save(std::string filename, std::function<double(double)> const& func_x);
    private:
        std::vector<std::vector<double>> InterpolateDistance(const WalkTrajectory& traj) {throw std::runtime_error("Not implemented");} //Interpolate trajectory to give time grid
        bool is_keep_trajectory() {return true;} //Do not store the full trajectory; just store the sample points

};

//Can deal with variable timestep
class WalkVariableTime : public Walk {
	public:
		WalkVariableTime(Lattice& latt0);
		//Save average of all trajectories ; a function convert distance to something is passed and averaged
		void save(std::string filename, std::function<double(double)> const& func_x);
		WalkTrajectory RunWalkOnce(double time_limit, double num_sample);
    private:
        std::vector<std::vector<double>> InterpolateDistance(const WalkTrajectory& traj); //Interpolate trajectory to give time grid
        bool is_keep_trajectory() {return false;} //Do not store the full trajectory; just store the sample points
        double EstimateSampleTime(double num_step_per_sample); //Estimate the temperature used to estimate the sample
};
