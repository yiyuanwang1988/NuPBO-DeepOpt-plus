#include "basis_pms.h"
#include "pms.h"
#include <sstream>

int main(int argc, char *argv[])
{
	// srand((unsigned)time(NULL));
	start_timing();
	Satlike s;
	stringstream ss;
	ss.str(argv[2]);
	int cutoff;
	ss >> cutoff;
	s.cutoff_time = cutoff;

	/*seed*/
	ss.clear();
	ss.str(argv[3]);
	int seed;
	ss >> seed;
	srand(seed);

	vector<int> init_solution;
	s.build_instance(argv[1]);

	s.local_search_with_decimation(init_solution, argv[1]);

	// s.simple_print();
	s.get_obj(argv[4]);
	s.print_best_solution(argv[1], argv[3]);
	s.free_memory();

	return 0;
}
