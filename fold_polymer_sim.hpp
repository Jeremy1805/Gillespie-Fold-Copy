#ifndef _PolySim_H
#define _PolySim_H

#include <vector>
#include <chrono>
#include "tsv_parser.hpp"

using namespace std;

template<class T, int N>
using boostmat = typename boost::multi_array<T,N>;

typedef boost::multi_array<int, 2> array_type;
typedef boost::multi_array_types::index_range range_t;

// Function to divide a vector/array by a scalar
template <typename Type>
vector<Type> div_scalar(vector<Type> &V, Type S);

// Function to perform a cumulative sum of the vector
template <typename Type>
vector<Type> cumsum(boostmat<Type,1> &V);

// Function to print arrays/vectors
template <class Type>
void print_arr(vector<Type> V);

// Timer: Specifically for Performance Measurements
class Timer
{
    private:
        chrono::high_resolution_clock::time_point start_time;
        chrono::high_resolution_clock::time_point stop_time;

    public:
        void start()
        {
            start_time = std::chrono::high_resolution_clock::now();
        }

        void stop()
        {
            stop_time = std::chrono::high_resolution_clock::now();
        }

        double measure()
        {
            using namespace std::chrono;
            return duration<float>
                (stop_time - start_time).count();
        }
};

// FoldDrivenSim: Class containing instructions on how to simulate
    // polymer copying coupled to polymer folding.

struct FoldDrivenSim {

    //**VARIABLES*//
    
    // Parameters
    double G_pol; 
    double k_copy;
    double k_fold;
    double g_fold;
    double g_u;
    double g_h;

    boostmat<double,2> energy_matrix{boost::extents[1][1]};
    int kinetic_flag;
    
    // Parameter-Derived
    boostmat<double,4> forward_rates{boost::extents[1][1][1][1]}; // 4-dim tensor with rate constants for polymer growth 
    boostmat<double,5> backward_rates{boost::extents[1][1][1][1][1]}; // 5-dim tensor with rate constants for monomer removal
     
    // States
    vector<int> fold;
    vector<int> copy;
    vector<int> template_s; // the s is only because 'template' is forbidden
    
    vector<int> complete_fold;
    vector<int> complete_copy;
    
    // State-Derived
    int polymer_size;
    int cpy_pos_start;

    //**FUNCTIONS**//

    // Constructor
    FoldDrivenSim();
    FoldDrivenSim(const FoldDrivenSim &tsv_parser);
    FoldDrivenSim& operator=(const FoldDrivenSim &tsv_parser);
    FoldDrivenSim(vector<int> initial_copy, vector<int> initial_template, vector<int> initial_fold, double G_pol_in, double k_copy_in,
        double k_fold_in,  double g_fold_in, double g_u_in, double g_h_in, boostmat<double,2> energy_matrix_in, int kinetic_flag_in); 
    FoldDrivenSim(TsvParser tsv_parser); 

    // Set Up
    int CheckInput(); // Checks input consistency and also returns the position of the first monomer to be added in the copy.
    void CreateKineticTensors(); // Create forward rates and backward rates objects

    // Simulation
    void Simulate(unsigned int seed);
    void DebugSimulator(unsigned int seed);
    void FoldCoilFraction(unsigned int seed, double t_lim, vector<int> &number_of_coil_states, vector<int> &tip_coil_states);
    void SpeedTest(unsigned int seed);
    vector<double> CreateWeightThresholds(); // Creating the initial probability thresholds
    
        
        // Simulation debug
        [[gnu::always_inline]] 
        inline void DebugNormalization( double partition_function, vector<double> weight_thresholds );
        [[gnu::always_inline]] 
        inline void DebugOutput( double rand_num, int newstate_flag, double partition_function, vector<double> &weight_thresholds,
            vector<int> &fold_run, vector<int> &copy_run, vector<int> &template_run, double microsecond_sleep_time );
        // Simulation subfunctions (short functions with potentially massive numbers of calls, inlined and arrays passed by reference for speed)
        [[gnu::always_inline]] 
        inline int GetContigHelix( int z, int a, int b, int d, int e, int f); // identify the size of the longest chain in the helix 
        [[gnu::always_inline]] 
        inline int SampleWeightThresholds( double rand_num, int weight_size, vector<double> &weight_thresholds);
        [[gnu::always_inline]] 
        inline int UpdateStateAndProb(int newstate_flag, int cpy_pos, int copy_alph_p_1, boostmat<double,5> &backward_rates_c, boostmat<double,4> &forward_rates_c,
            double param_Kfold, double param_Gfold, double param_Gu, double param_Gh, vector<int>& fold_r, vector<int>& copy_r, vector<int> &template_r, vector<double> &weight_thresholds, double &partition_func);        
        [[gnu::always_inline]] 
        inline int UpdateHbondDebug(int newstate_flag, int cpy_pos, int copy_alph_p_1, boostmat<double,5> &backward_rates_c, boostmat<double,4> &forward_rates_c,
            double param_Kfold, double param_Gfold, double param_Gu, double param_Gh, vector<int>& fold_r, vector<int>& copy_r, vector<int> &template_r, vector<double> &weight_thresholds, double &partition_func, int &hydrogen_count);

    // Utility
    void PrintKineticTensors(); // Print forward rates and backward rates
    int EvalHydrobonds(vector<int> fold); // Hydrogen bond counting for debug only
    inline int EvalCoil( vector<int> fold );
    double GetAccuracy();
};

#endif
