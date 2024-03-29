Kinetic Monte Carlo Random Walk Simulation
==========================================


Description
------------------------------------

This code computes a random walk simulation of polaron hopping with kinetic monte carlo sampling as described in this article:

[F. Wu, and Y. Ping, *Journal of Material Chemistry A* **6**, 20025-20036 (2018).](https://pubs.rsc.org/en/content/articlehtml/2018/ta/c8ta07437b "Combining Landau–Zener theory and kinetic Monte Carlo sampling for small polaron mobility of doped BiVO4 from first-principles")


Run the program
------------------------------------

`randomwalk` is the main program.

Two files required, $PREFIX.in and $PREFIX.latt. 

Run `randomwalk  $PREFIX` to run the program.

`calc_activation.py` is the script to automate the calculation and post-processing.

* To verify the results, with `-n $N` option, this program will run the simulation for $N times, and compare results to get error estimation. (Note this is individual to `num_trajectory` in the input file)
    
* This script supports parallel runs on a single node with parameter "-np" (suggest `np = n`)

`kairay.sh` is an example sbatch script on Kairay to run 3200 and 12800 iterations to show the convergence.


Input file description
------------------------------------

*.in :  consists of only numbers

num_step_per_sample :  run simulation for num_step_per_sample time steps for each L(t)-t data

num_sample  : record L(t)-t data in one simulation for num_sample times

num_trajectory : repeat the simulation of num_trajectory times and take the average

number of temperatures

temperature 1

temperature 2

...


Known issues
------------------------------------

The time step is estimated from one test run; so if there are sites with very different barriers, the estimated time would be either too much for fast hopping or too less for slow hopping. If the first happens then the program will go out of memory. 

Author(s)
------------------------------------
Feng Wu

Tyler Smart

Stefano Falletta 

