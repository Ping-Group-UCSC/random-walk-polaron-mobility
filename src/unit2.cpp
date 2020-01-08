/*Test some functions
 */
#include "param.h"
#include "walk.h"

/*Test random functions
 */

int main(int argc, char* argv[]){
	std::fstream fs("bivo4.latt", std::fstream::in);
	if (!fs.good()) {
		throw std::runtime_error("Cannot open input file : latt");
	}
	LatticeGT latt = LatticeGT();
	fs >> latt;
	fs.close();
    WalkVariableTime walk = WalkVariableTime(latt);
    walk.random.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    latt.SetTemperature(300);
    std::map<int, int> grid_visit;
    grid_visit.clear();

    for (int i = 0 ; i < 100000 ; i ++){
//       std::cout << walk.RandomSelectPath(latt.grids[0]).index_grid << " ";
         int ig = walk.RandomSelectPath(latt.grids[0]).index_grid;
         if (grid_visit.find(ig) == grid_visit.end()){
             grid_visit[ig] = 0;
         }
         grid_visit[ig] += 1;
    }
    for (auto const & path : latt.grids[0].paths){
        std::cout << path.k << " " << latt.grids[path.index_grid].name << " " << grid_visit[path.index_grid] << std::endl;
    }

    std::cout << std::endl;
//  int n = 3;
//  for (int i = 0 ; i < 10000000000 ; i++){
//      if (walk.RandomSelect(n) == n){
//          std::cout << "Meet random number == maximum" << std::endl;
//      }
//  }
}
