#include "boost/multi_array.hpp"

#include "fold_polymer_sim.hpp"
#include "tsv_parser.hpp"

template<class T, int N>
using boostmat = typename boost::multi_array<T,N>;

typedef boost::multi_array<int, 2> array_type;
typedef boost::multi_array_types::index_range range_t;

using namespace std;


struct BatchSimulator {
    TsvParser tsv_parser;
    FoldDrivenSim fold_driven_sim;
    
    BatchSimulator(string);
    void SimulateAll();
};


