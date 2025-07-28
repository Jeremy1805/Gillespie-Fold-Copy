#include <iostream>
#include <fstream>
#include <vector>

#include "boost/multi_array.hpp"
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"

#include "batch_simulator.hpp"
#include "tsv_parser.cpp"
#include "fold_polymer_sim.cpp"

template<class T, int N>
using boostmat = typename boost::multi_array<T,N>;

typedef boost::multi_array<int, 2> array_type;
typedef boost::multi_array_types::index_range range_t;

using namespace std;

// Function for converting a 2-dim boost matrix into
// string "a,b;c,d" where , deliniates columns
// and ; delineates rows

// Remove the first row and first column corresponding to zeros

BatchSimulator::BatchSimulator(string input_string) {
    tsv_parser = TsvParser(input_string);
    fold_driven_sim = FoldDrivenSim(tsv_parser);
}

void BatchSimulator::SimulateAll() {

    vector< vector<int> > copies;
    vector< vector<int> > folds;
    vector<double> accuracy_ls;
    double mean_acc;

    double sim_num = double(tsv_parser.simulation_repeat);

    for ( unsigned int seed : tsv_parser.random_seed) {
        cout << seed << endl;
        fold_driven_sim.Simulate(seed);

        copies.push_back(fold_driven_sim.complete_copy);
        folds.push_back(fold_driven_sim.complete_fold);
        double acc = fold_driven_sim.GetAccuracy();
        accuracy_ls.push_back( acc );
        mean_acc = mean_acc + acc/sim_num;
    } 

    tsv_parser.output_copy_list = copies;
    tsv_parser.output_folds_list = folds;
    tsv_parser.accuracy_list = accuracy_ls;
    tsv_parser.mean_accuracy = mean_acc;
}
