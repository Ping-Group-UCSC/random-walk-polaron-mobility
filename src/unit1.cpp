/*Test some functions
 */
#include "param.h"
#include "walk.h"


int main(int argc, char* argv[]){
	std::fstream fs("bivo4-16.latt", std::fstream::in);
	if (!fs.good()) {
		throw std::runtime_error("Cannot open input file : latt");
	}
	LatticeGT latt = LatticeGT();
	fs >> latt;
	fs.close();
    WalkVariableTime walk = WalkVariableTime(latt);
    walk.random.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    for (int i = 0 ; i < 100 ; i ++){
        std::cout << walk.RandomSelect(3) << " ";
    }
    std::cout << std::endl;
    int n = 3;
    for (int i = 0 ; i < 10000000000 ; i++){
        if (walk.RandomSelect(n) == n){
            std::cout << "Meet random number == maximum" << std::endl;
        }
    }
}
