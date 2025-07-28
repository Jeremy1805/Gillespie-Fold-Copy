#ifndef _PolySim
#define _PolySim

// In this version of the program, we keep track of the normalization constant separate from the underlying partition function. 
// However, due to the use of exponentials, this results in a drift in the normalization factor away from the vector. 

#include <iostream>
#include <numeric>
#include <iterator>
#include <vector>
#include <unordered_map>
#include <math.h>
#include <random>
#include <chrono>
#include<unistd.h>
#include <fstream>

#include "boost/multi_array.hpp"
#include "tsv_parser.cpp"
#include "fold_polymer_sim.hpp"

using namespace std;

// cpp file contains definitions for the abstract GillespieSim class, but
// not the specific implementations 

/******************
 GENERIC FUNCTIONS
 ******************/

// Function to divide a vector/array by a scalar
template <typename Type>
vector<Type> div_scalar(vector<Type> &V, Type S) {
    vector<Type> Ans;
    Ans.reserve(V.size());
 
    for(typename vector<Type>::iterator it = V.begin(); it != V.end(); ++it) {
        Ans.push_back(*it/S);
    }
    return Ans;
}

// Function to perform a cumulative sum of the vector
template <typename Type>
vector<Type> cumsum(boostmat<Type,1> &V) {
    vector<Type> Ans;
    
    Type running_total = 0;
    for(int i = 0; i != V.shape()[0]; ++i) {
        running_total += V[i];
        Ans.push_back(running_total);
    }

    return Ans;
}

// Function to print arrays/vectors
template <class Type>
void print_arr(vector<Type> V) {
    for(typename vector<Type>::iterator it = V.begin(); it != V.end(); ++it) {
        if (it != V.end()-1){
            cout << *it << " ";
	} else {
	    cout << *it;
	}
    }
    cout << flush;
}

string stringify(vector<int> V) {
    string out="";
    for(typename vector<int>::iterator it = V.begin(); it != V.end(); ++it) {
        out.append(to_string(*it));            
	}
    return(out);
}

// Function to print arrays/vectors
template <class Type>
void print_arr_with_idx(vector<Type> V, int idx) {
    for(typename vector<Type>::iterator it = V.begin(); it != V.begin()+idx; ++it) {
        if (it != V.end()-1){
            cout << *it << " ";
	} else {
	    cout << *it;
	}
    }
    cout << flush;
}

/******************
 Class: FoldDrivenSim
 ******************/

// Abstract class containing instructions on how to simulate
// as well as converting from a tensor of reaction rates
// to a tensor of probabilities and accessible states

//**CONSTRUCTOR**//
FoldDrivenSim::FoldDrivenSim(){};
FoldDrivenSim::FoldDrivenSim( vector<int> initial_copy, vector<int> initial_template, vector<int> initial_fold, double G_pol_in, double k_copy_in,
    double k_fold_in, double g_fold_in, double g_u_in, double g_h_in, boostmat<double,2> energy_matrix_in, int kinetic_flag_in) {
    
    // Parameters
    G_pol = G_pol_in;
    k_copy = k_copy_in;
    k_fold = k_fold_in;
    g_fold = g_fold_in;
    g_u = g_u_in;
    g_h = g_h_in;
    energy_matrix.resize(boost::extents[energy_matrix_in.shape()[0]][energy_matrix_in.shape()[1]]);
    energy_matrix = energy_matrix_in;
    kinetic_flag = kinetic_flag_in;
   
    // Simulation State Set Up
    fold = initial_fold;
    copy = initial_copy;
    template_s = initial_template;
    cpy_pos_start = CheckInput();

    // Obtain Sizes 
    polymer_size = template_s.size();
    CreateKineticTensors();   
    PrintKineticTensors();
    cout << endl;
}

FoldDrivenSim::FoldDrivenSim( TsvParser tsv_parser) {

    // Parameters
    G_pol = tsv_parser.G_pol;
    k_copy = tsv_parser.k_copy;
    k_fold = tsv_parser.k_fold;
    g_fold = tsv_parser.G_fold;
    g_u = tsv_parser.G_u;
    g_h = tsv_parser.G_h;
    energy_matrix.resize(boost::extents[tsv_parser.energy_matrix.shape()[0]][tsv_parser.energy_matrix.shape()[1]]);
    energy_matrix = tsv_parser.energy_matrix;
    kinetic_flag = tsv_parser.kinetic_flag;
   
    // Simulation State Set Up
    fold = tsv_parser.initial_fold;
    copy = tsv_parser.initial_copy;
    template_s = tsv_parser.initial_template;
    cpy_pos_start = CheckInput();

    // Obtain Sizes 
    polymer_size = template_s.size();
    CreateKineticTensors();   
    PrintKineticTensors();
    cout << endl;
}

FoldDrivenSim::FoldDrivenSim(const FoldDrivenSim &tsv_parser) {
    G_pol = tsv_parser.G_pol;
    k_copy = tsv_parser.k_copy;
    k_fold = tsv_parser.k_fold;
    g_fold = tsv_parser.g_fold;
    g_u = tsv_parser.g_u;
    g_h = tsv_parser.g_h;
    energy_matrix.resize(boost::extents[tsv_parser.energy_matrix.shape()[0]][tsv_parser.energy_matrix.shape()[1]]);
    energy_matrix = tsv_parser.energy_matrix;
    kinetic_flag = tsv_parser.kinetic_flag;
   
    // Simulation State Set Up
    fold = tsv_parser.fold;
    copy = tsv_parser.copy;
    template_s = tsv_parser.template_s;
    cpy_pos_start = CheckInput();

    // Obtain Sizes 
    polymer_size = template_s.size();
    CreateKineticTensors(); 
}

FoldDrivenSim& FoldDrivenSim::operator=(const FoldDrivenSim &tsv_parser) {
    G_pol = tsv_parser.G_pol;
    k_copy = tsv_parser.k_copy;
    k_fold = tsv_parser.k_fold;
    g_fold = tsv_parser.g_fold;
    g_u = tsv_parser.g_u;
    g_h = tsv_parser.g_h;
    energy_matrix.resize(boost::extents[tsv_parser.energy_matrix.shape()[0]][tsv_parser.energy_matrix.shape()[1]]);
    energy_matrix = tsv_parser.energy_matrix;
    kinetic_flag = tsv_parser.kinetic_flag;
   
    // Simulation State Set Up
    fold = tsv_parser.fold;
    copy = tsv_parser.copy;
    template_s = tsv_parser.template_s;
    cpy_pos_start = CheckInput();

    // Obtain Sizes 
    polymer_size = template_s.size();
    CreateKineticTensors(); 
    return *this;
}
//**SET UP**//

// Function to Check the Input Consistency of the Polymer Copy, Template, Fold and Energy Matrix. 
    // Also returns the position of the first monomer to be added in the copy. 
int FoldDrivenSim::CheckInput(){
        
    // First check that the sizes of the vectors are consistent
    if ( fold.size() != copy.size() || copy.size() != template_s.size() ) {
        throw runtime_error("ERROR: INCONSISTENT VECTOR SIZES!");
    }

    // Next, check that the very first character is always a 0 in all cases 
    if ( fold[0] != 0 || copy[0] != 0 || template_s[0] != 0 ) {
        throw runtime_error("ERROR: INCORRECT STARTING CHARACTERS!");
    }

    // Triple padding makes checks for fold flipping easier, so check that the next two characters are always a 0 in all cases 
    if ( fold[1] != 0 || copy[1] != 0 || template_s[1] != 0 ) {
        throw runtime_error("ERROR: INCORRECT STARTING CHARACTERS!");
    }

    if ( fold[2] != 0 || copy[2] != 0 || template_s[2] != 0 ) {
        throw runtime_error("ERROR: INCORRECT STARTING CHARACTERS!");
    }
    
    // Next, check that the very last character is always a 0 in all cases 
    if ( fold[fold.size() - 1] != 0 || copy[ copy.size() - 1 ] != 0 || template_s[ template_s.size() - 1 ] != 0 ) {
        throw runtime_error("ERROR: INCORRECT ENDING CHARACTERS!");
    }

    // Triple padding again 
    if ( fold[fold.size() - 2] != 0 || copy[ copy.size() - 2 ] != 0 || template_s[ template_s.size() - 2 ] != 0 ) {
        throw runtime_error("ERROR: INCORRECT ENDING CHARACTERS!");
    }

    if ( fold[fold.size() - 3] != 0 || copy[ copy.size() - 3 ] != 0 || template_s[ template_s.size() - 3 ] != 0 ) {
        throw runtime_error("ERROR: INCORRECT ENDING CHARACTERS!");
    }

    // Identify ranges that the alphabet can take based on the provided energy matrix
    int largest_copy = energy_matrix.shape()[0] - 1; 
    int largest_template = energy_matrix.shape()[1] - 1;
    
    // Then, we check that the last non-zero index(besides location 0) is equal between the copy and fold 
        // and check that all characters are within bounds

    int first_zero_idx;
    for (int cpy_idx = 3; cpy_idx < copy.size(); cpy_idx++) {
        
        if ( copy[cpy_idx] < 0 || copy[cpy_idx] > largest_copy ) {
            throw runtime_error("ERROR: OUT OF BOUND MONOMER IN THE COPY!");
        }

        if ( copy[cpy_idx] == 0 ) {
            first_zero_idx = cpy_idx;
            break;
        }
    }

    for (int cpy_idx = first_zero_idx; cpy_idx < copy.size(); cpy_idx++) {
        if ( copy[cpy_idx] !=0 ) {
            throw runtime_error("ERROR: DISCONTIGUOUS COPY!");
        }
    }
    
    for (int fld_idx = 3; fld_idx < first_zero_idx; fld_idx++) {
        if ( fold[fld_idx] < 0 || fold[fld_idx] > 2 ) {
            throw runtime_error("ERROR: OUT OF BOUND MONOMER IN THE FOLD!");
        }
        
        if ( fold[fld_idx] == 0 ) {
            throw runtime_error("ERROR: INCONSISTENT COPY AND FOLD STATES!");
        }
    }

    for (int fld_idx = first_zero_idx; fld_idx < fold.size(); fld_idx++) {
        if ( fold[fld_idx] !=0 ) {
            throw runtime_error("ERROR: DISCONTIGUOUS FOLD!");
        }
    }

    for (int tpt_idx = 3; tpt_idx < template_s.size() - 3; tpt_idx++) {
        if ( template_s[tpt_idx] <= 0 || template_s[tpt_idx] > largest_template ) {
            throw runtime_error("ERROR: OUT OF BOUND MONOMER IN THE TEMPLATE!");
        }
    }
    
    if ((fold[first_zero_idx - 1] != 2) && (first_zero_idx > 3)) {
        throw runtime_error("ERROR: Tip monomer must be in 'coil' state (not a true coil, but a representation of the monomer bound to template state )");
    }

    if ((fold[3] != 2)) {
        throw runtime_error("ERROR:  There must be a start monomer, and start monomer must be in 'coil' state (not a real monomer as it cannot flip/be removed, just a boundaryk)");
    }


    return(first_zero_idx - 1);
}

void FoldDrivenSim::CreateKineticTensors(){
    
    int copy_alphabet_size = energy_matrix.shape()[0];
    int template_alphabet_size = energy_matrix.shape()[1];
    
    // Note previous and next have different interpretations for forward and backward rates
        // In the forward case, next is the monomer to be added and previous is the position before it
        // In the backward case, next is the monomer to be removed and previous is the position before it

    // Dimensions are:
        // 0. The copy monomer in the previous position
        // 1. The template monomer in the previous position
        // 2. The copy monomer in the next position
        // 3. The template monomer in the next position
        // 4. The state of the fold in the next position (only relevant in the backward rate calculation)

    forward_rates.resize(boost::extents[copy_alphabet_size][template_alphabet_size][copy_alphabet_size][template_alphabet_size]);
    backward_rates.resize(boost::extents[copy_alphabet_size][template_alphabet_size][copy_alphabet_size][template_alphabet_size][3]);
    
    double forward_rate;
    double backward_rate;
    
    // Note: Index starts at 1 as 0 is reserved for the 'empty' location
    for (int copy_back = 1; copy_back != copy_alphabet_size; copy_back++){
        for (int template_back = 1; template_back != template_alphabet_size; template_back++){
            for (int copy_next = 1; copy_next != copy_alphabet_size; copy_next++){
                for (int template_next = 1; template_next != template_alphabet_size; template_next++){
                    if ( kinetic_flag == 0 ){
                        // Homogeneous case uses energies from energy matrix to determine 
                        // forward rates
                        forward_rate = k_copy*exp(energy_matrix[copy_next][template_next]);
                    } else if ( kinetic_flag == 1 ){
                        // Inhomogeneous case always sets forward rate to k_copy
                        forward_rate = k_copy;
                    } else {
                        throw runtime_error("ERROR: UNRECOGNIZED KINETIC FLAG!");
                    }

                    backward_rate = exp(-G_pol
                            - energy_matrix[copy_next][template_next]
                        + energy_matrix[copy_back][template_back])*forward_rate;

                        forward_rates[copy_back][template_back][copy_next][template_next] = forward_rate;    
                        backward_rates[copy_back][template_back][copy_next][template_next][2] = backward_rate;
                    }
                }	    
            }
    }
}

//*SIMULATION*//
void FoldDrivenSim::Simulate(unsigned int seed){
    mt19937_64 mersenne_generator(seed);

    // initialize the weights and partition function
    vector<double> weight_thresholds = CreateWeightThresholds();

    double partition_function = std::accumulate(weight_thresholds.begin(), weight_thresholds.end(),
                                    decltype(weight_thresholds)::value_type(0));

    // For optimization, we keep a local copy of states and parameters
    vector<int> fold_run = fold;
    vector<int> copy_run = copy;
    vector<int> template_run = template_s;
    
    boostmat<double,5> backward_rates_copy;
    backward_rates_copy.resize(boost::extents[ backward_rates.shape()[0] ][ backward_rates.shape()[1] ][ backward_rates.shape()[2] ][ backward_rates.shape()[3] ][ backward_rates.shape()[4] ]);
    backward_rates_copy = backward_rates;
    boostmat<double,4> forward_rates_copy;
    forward_rates_copy.resize(boost::extents[ forward_rates.shape()[0] ][ forward_rates.shape()[1] ][ forward_rates.shape()[2] ][ forward_rates.shape()[3] ]);
    forward_rates_copy = forward_rates;

    int final_pos = polymer_size - 4;
    int cpy_pos = cpy_pos_start;
    int polymer_sz = polymer_size;
    int copy_alph_p_1 = energy_matrix.shape()[0] + 1;
    int weight_size = weight_thresholds.size();

    double param_Kfold = k_fold;
    double param_Gfold = g_fold;
    double param_Gu = g_u;
    double param_Gh = g_h;

    double rand_num;
    int newstate_flag;
    while (cpy_pos < final_pos){   
        rand_num = uniform_real_distribution<double>(0.0,partition_function)(mersenne_generator); 

        newstate_flag = SampleWeightThresholds(rand_num,weight_size,weight_thresholds);
        // Get the index of the next state in the weight vector using a naive gillespie method
        // See the comments on the function CreateWrightThresholds to see how state indices are partitioned

        cpy_pos = UpdateStateAndProb(newstate_flag, cpy_pos, copy_alph_p_1, backward_rates_copy, forward_rates_copy,
            param_Kfold, param_Gfold, param_Gu, param_Gh, fold_run, copy_run, template_run, weight_thresholds, partition_function);
        // Update state and Weight thresholds. Return the new position of the copy (no change if the state change is a change in the fold state instead)

        // Debug only: check normalization
            // DebugNormalization(partition_function, weight_thresholds);
        /* Debug only: Print
            DebugOutput( rand_num, newstate_flag, partition_function, weight_thresholds,
                fold_run, copy_run, template_run, 10000 );
        */
    }
    complete_copy = copy_run;
    complete_fold = fold_run;
}

void FoldDrivenSim::DebugSimulator(unsigned int seed){
    mt19937_64 mersenne_generator(seed);
    
    // initialize the weights and partition function
    vector<double> weight_thresholds = CreateWeightThresholds();
    double partition_function = std::accumulate(weight_thresholds.begin(), weight_thresholds.end(),
                                    decltype(weight_thresholds)::value_type(0));

    // For optimization, we keep a local copy of states and parameters
    vector<int> fold_run = fold;
    vector<int> copy_run = copy;
    vector<int> template_run = template_s;

    boostmat<double,5> backward_rates_copy;
    backward_rates_copy.resize(boost::extents[ backward_rates.shape()[0] ][ backward_rates.shape()[1] ][ backward_rates.shape()[2] ][ backward_rates.shape()[3] ][ backward_rates.shape()[4] ]);
    backward_rates_copy = backward_rates;
    boostmat<double,4> forward_rates_copy;
    forward_rates_copy.resize(boost::extents[ forward_rates.shape()[0] ][ forward_rates.shape()[1] ][ forward_rates.shape()[2] ][ forward_rates.shape()[3] ]);
    forward_rates_copy = forward_rates;

    int final_pos = polymer_size - 4;
    int cpy_pos = cpy_pos_start;
    int polymer_sz = polymer_size;
    int copy_alph_p_1 = energy_matrix.shape()[0] + 1;
    int weight_size = weight_thresholds.size();

    double param_Kfold = k_fold;
    double param_Gfold = g_fold;
    double param_Gu = g_u;
    double param_Gh = g_h;

    double rand_num;
    int newstate_flag;
    int hydrogen_count = 0;

    //ofstream fold_log;
    //fold_log.open ("fold_log.txt");
    while (cpy_pos < final_pos){   
        rand_num = uniform_real_distribution<double>(0.0,partition_function)(mersenne_generator); 

        newstate_flag = SampleWeightThresholds(rand_num,weight_size,weight_thresholds);
        // Get the index of the next state in the weight vector using a naive gillespie method
        // See the comments on the function CreateWrightThresholds to see how state indices are partitioned

        cpy_pos = UpdateHbondDebug(newstate_flag, cpy_pos, copy_alph_p_1, backward_rates_copy, forward_rates_copy,
            param_Kfold, param_Gfold, param_Gu, param_Gh, fold_run, copy_run, template_run, weight_thresholds, partition_function, hydrogen_count);
        // Update state and Weight thresholds. Return the new position of the copy (no change if the state change is a change in the fold state instead)
        
        //fold_log << stringify(fold_run);

        //Debug only: check normalization
            DebugNormalization(partition_function, weight_thresholds);
        // Debug only: Print
            DebugOutput( rand_num, newstate_flag, partition_function, weight_thresholds,
                fold_run, copy_run, template_run, 10000 );
            
            cout << "Hydrogen Count from Simulation: " << hydrogen_count << endl;
            cout << "Hydrogen Count from Evaluation: " << EvalHydrobonds(fold_run) << endl;

            cout << endl;
            /*
            if (hydrogen_count != EvalHydrobonds(fold_run)) {
                throw runtime_error("Hydrogen bond counts not equal!");
            }*/
    }
    complete_copy = copy_run;
    complete_fold = fold_run;
    //fold_log.close();
}

// Simulation of the fold only, all copy related rates set to 0
    // The end goal is to obtain a vector of the number of 'coil' states at each point in time
    // as well as the state of the tip at each point in time
void FoldDrivenSim::FoldCoilFraction(unsigned int seed, double t_lim, vector<int> &out_coil_number, vector<int> &out_tip_states){
    mt19937_64 mersenne_generator(seed);
    
    if ( k_copy != 0 ) {
        throw runtime_error("Fold Simulation requested, please set k_copy to 0!");
    }

    // initialize the weights and partition function
    vector<double> weight_thresholds = CreateWeightThresholds();
    double partition_function = std::accumulate(weight_thresholds.begin(), weight_thresholds.end(),
                                    decltype(weight_thresholds)::value_type(0));

    // For optimization, we keep a local copy of states and parameters
    vector<int> fold_run = fold;
    vector<int> copy_run = copy;
    vector<int> template_run = template_s;

    boostmat<double,5> backward_rates_copy;
    backward_rates_copy.resize(boost::extents[ backward_rates.shape()[0] ][ backward_rates.shape()[1] ][ backward_rates.shape()[2] ][ backward_rates.shape()[3] ][ backward_rates.shape()[4] ]);
    // backward_rates_copy = backward_rates;
        // all copy rates set to 0
    boostmat<double,4> forward_rates_copy;
    forward_rates_copy.resize(boost::extents[ forward_rates.shape()[0] ][ forward_rates.shape()[1] ][ forward_rates.shape()[2] ][ forward_rates.shape()[3] ]);
    // forward_rates_copy = forward_rates;

    int final_pos = polymer_size - 4;
    int cpy_pos = cpy_pos_start;
    int polymer_sz = polymer_size;
    int copy_alph_p_1 = energy_matrix.shape()[0] + 1;
    int weight_size = weight_thresholds.size();

    double param_Kfold = k_fold;
    double param_Gfold = g_fold;
    double param_Gu = g_u;
    double param_Gh = g_h;

    double rand_num;
    double rand_num_2;

    int newstate_flag;
    
    vector<int> coil_number;
    vector<int> tip_states;

    double t = 0;
    double dt = 0;
    // reserve a number of vector elements with our best guess of the number of transitions 
    coil_number.reserve( int(t_lim/partition_function) );
    tip_states.reserve( int(t_lim/partition_function) );

    int running_coil_number = 0;

    for (int m: fold_run) {
        if (m == 2) {
            running_coil_number++;
        }
    }
    
    //ofstream fold_log;
    //fold_log.open ("fold_log.txt");
    // when sampling for distributions, it is vital that we stop at a specific time rather than a specific number of steps
    while (t < t_lim){   
        
        coil_number.push_back( running_coil_number );
        tip_states.push_back( fold_run[ cpy_pos - 1 ] );
        
        rand_num = uniform_real_distribution<double>(0.0,partition_function)(mersenne_generator); 
        rand_num_2 = uniform_real_distribution<double>(0.0,1.0)(mersenne_generator);
        
        dt = - log( 1 - rand_num_2 )/partition_function;

        t = t + dt;

        newstate_flag = SampleWeightThresholds(rand_num,weight_size,weight_thresholds);
        
        running_coil_number = EvalCoil(fold_run);
        
        cpy_pos = UpdateStateAndProb(newstate_flag, cpy_pos, copy_alph_p_1, backward_rates_copy, forward_rates_copy,
            param_Kfold, param_Gfold, param_Gu, param_Gh, fold_run, copy_run, template_run, weight_thresholds, partition_function);
        
        //fold_log << stringify(fold_run) << ";" << dt << "\n";

        // Debug only: check normalization
            // DebugNormalization(partition_function, weight_thresholds);
        /* Debug only: Print
            DebugOutput( rand_num, newstate_flag, partition_function, weight_thresholds,
                fold_run, copy_run, template_run, 10000 );
        */
    }
    //fold_log.close();
    
    out_coil_number = coil_number;
    out_tip_states = tip_states;
}

//*SIMULATION*//
void FoldDrivenSim::SpeedTest(unsigned int seed){
    mt19937_64 mersenne_generator(seed);
    
    // initialize the weights and partition function
    vector<double> weight_thresholds = CreateWeightThresholds();
    double partition_function = std::accumulate(weight_thresholds.begin(), weight_thresholds.end(),
                                    decltype(weight_thresholds)::value_type(0));

    // For optimization, we keep a local copy of states and parameters
    vector<int> fold_run = fold;
    vector<int> copy_run = copy;
    vector<int> template_run = template_s;

    boostmat<double,5> backward_rates_copy;
    backward_rates_copy.resize(boost::extents[ backward_rates.shape()[0] ][ backward_rates.shape()[1] ][ backward_rates.shape()[2] ][ backward_rates.shape()[3] ][ backward_rates.shape()[4] ]);
    backward_rates_copy = backward_rates;
    boostmat<double,4> forward_rates_copy;
    forward_rates_copy.resize(boost::extents[ forward_rates.shape()[0] ][ forward_rates.shape()[1] ][ forward_rates.shape()[2] ][ forward_rates.shape()[3] ]);
    forward_rates_copy = forward_rates;

    int final_pos = polymer_size - 4;
    int cpy_pos = cpy_pos_start;
    int polymer_sz = polymer_size;
    int copy_alph_p_1 = energy_matrix.shape()[0] + 1;
    int weight_size = weight_thresholds.size();

    double param_Kfold = k_fold;
    double param_Gfold = g_fold;
    double param_Gu = g_u;
    double param_Gh = g_h;

    double rand_num;
    int newstate_flag;

    bool dummy;

    Timer t; // this is not simulation time, but a class for measuring computation time
    t.start();

    for(int k = 0; k != 2000000; k++){   
        rand_num = uniform_real_distribution<double>(0.0,partition_function)(mersenne_generator); 

        newstate_flag = SampleWeightThresholds(rand_num,weight_size,weight_thresholds);
        // Determine the current state (defined by the 5 numbers)
            // 1. copy_back 2.template_back 3.copy_current 4. template_current 5. template_next

        cpy_pos = UpdateStateAndProb(newstate_flag, cpy_pos, copy_alph_p_1, backward_rates_copy, forward_rates_copy,
            param_Kfold, param_Gfold, param_Gu, param_Gh, fold_run, copy_run, template_run, weight_thresholds, partition_function);
        
        dummy = cpy_pos < final_pos;
    }

    t.stop();
    double measured = t.measure();
    
    cout << measured << endl;
}

vector<double> FoldDrivenSim::CreateWeightThresholds(){
    // Entry 0 of the vector is the weight of monomer removal
    // Entries 1 to C, where C is the number of possible copy monomers, are reservedfor monomer addition
    // Entries C+1 to C+N, where N is the polymer length, are reserved for fold flipping 
    
    int copy_alph_size = energy_matrix.shape()[0];

    vector<double> output_weights(polymer_size+copy_alph_size+1,0);
    for ( int i = 0; i < copy_alph_size; i++ ) {
        output_weights[1+i] = forward_rates[ copy[cpy_pos_start] ][ template_s[cpy_pos_start] ]
            [ i ][ template_s[cpy_pos_start + 1] ];
    }

    for ( int i = 3; i < polymer_size-2; i++ ) {
        if (i == 3) {
            output_weights[1 + copy_alph_size + i] = 0;
        } else if (i == cpy_pos_start) {
            output_weights[1 + copy_alph_size + i] = 0; // The copy position cannot flip
        } else if (fold[i] == 1){
            output_weights[1 + copy_alph_size + i] = k_fold*exp( -GetContigHelix( fold[i-3], fold[i-2], fold[i-1], fold[i+1], fold[i+2], fold[i+3])*(g_h-g_u) - g_fold  );
        } else if (fold[i] == 2) {
            output_weights[1 + copy_alph_size + i] = k_fold*exp( GetContigHelix( fold[i-3],fold[i-2], fold[i-1], fold[i+1], fold[i+2], fold[i+3]  )*g_u );
        } else {
            output_weights[1 + copy_alph_size + i] = 0;
        }
    }
    
    return(output_weights);
    
}

inline void FoldDrivenSim::DebugNormalization( double partition_function, vector<double> weight_thresholds ){
    double truesum = std::accumulate(weight_thresholds.begin(), weight_thresholds.end(),
            decltype(weight_thresholds)::value_type(0));
            
    if (truesum - partition_function > 1e-12) {
        cout << "Rounding Error: " << truesum - partition_function << endl;
        throw runtime_error("Normalization factor and weight vector out of sync!");
    }
}

inline void FoldDrivenSim::DebugOutput( double rand_num, int newstate_flag, double partition_function, vector<double> &weight_thresholds,
        vector<int> &fold_run, vector<int> &copy_run, vector<int> &template_run, double microsecond_sleep_time ){
    cout << "random: " << rand_num << endl;
    cout << "newstate_flag: " << newstate_flag << endl;
    cout << endl;

    print_arr(weight_thresholds);
    cout << endl;
    cout << "normalization factor: " << partition_function << endl; 
    cout << "SIMULATION: " << endl;
    print_arr(fold_run); cout << endl;
    print_arr(copy_run); cout << endl;
    print_arr(template_run); cout << endl;
                       
    usleep(microsecond_sleep_time);
}

inline int FoldDrivenSim::GetContigHelix(int z, int a, int b, int d, int e, int f){
    // Note function automatically assumes c is equal to 1
    int bonds = 0;
    if (z > 0 && d > 0) {
        bonds = bonds + int( (a == 1) && (b == 1) );
    }
    if ( a > 0 && e > 0) {
        bonds = bonds + int( (b == 1) && (d == 1) );
    }
    if ( b > 0 && f > 0) {
        bonds = bonds + int( (d == 1) && (e == 1) );
    }
    return(bonds);
}

inline int FoldDrivenSim::SampleWeightThresholds( double rand_num, int weight_size, vector<double> &weight_thresholds){
    double s = 0;
    for ( int i = 0; i < weight_size; i++ ) {
        s = s + weight_thresholds[i];
        if (rand_num < s) {
            return(i);
        }
    }
    return(weight_size-1); // This is possible as a result of rounding errors, 
        //but the magnitude of these errors is expected to be no greater than 1e-12 and hence negligible
}

inline int FoldDrivenSim::UpdateStateAndProb(int newstate_flag, int cpy_pos, int copy_alph_p_1, boostmat<double,5> &backward_rates_c, boostmat<double,4> &forward_rates_c,
            double param_Kfold, double param_Gfold, double param_Gu, double param_Gh, vector<int> &fold_r, vector<int> &copy_r, vector<int> &template_r, 
            vector<double> &weight_thresholds, double &partition_func){
    
    double diff;

    if (newstate_flag == 0) {
        // A monomer is removed.

        // STATE CHANGES: remove folds and copy 
        fold_r[cpy_pos] = 0;
        copy_r[cpy_pos] = 0;

        // The monomer behind the removed monomer can no longer flip as it is now bound to the copy, set the weight to zero.
        // Update the normalization factor first
        partition_func = partition_func - weight_thresholds[copy_alph_p_1+cpy_pos-1]; //each time a weight is updated, 
            // the normalization factor is also updated. The idea is to avoid having to sum up the whole vector each time
        weight_thresholds[copy_alph_p_1+cpy_pos-1] = 0;

        int new_cpy_pos = cpy_pos-1;
       
        // update monomer removal weights
        diff = backward_rates_c[ copy_r[new_cpy_pos-1] ][ template_r[new_cpy_pos-1] ]
            [ copy_r[new_cpy_pos] ][ template_r[ new_cpy_pos ] ][ fold_r[new_cpy_pos-1]] - weight_thresholds[0];
        
        weight_thresholds[0] = weight_thresholds[0] + diff;
        partition_func = partition_func + diff;

        for ( int i = 0; i < copy_alph_p_1 - 1; i++ ) {
            // update monomer addition weights
            diff = forward_rates_c[ copy_r[new_cpy_pos] ][ template_r[new_cpy_pos] ]
                [ i ][ template_r[new_cpy_pos + 1] ] - weight_thresholds[i+1];

            weight_thresholds[i+1] = weight_thresholds[i+1] + diff;
            partition_func = partition_func + diff;
        };
        
        return(new_cpy_pos);
    } else if (newstate_flag < copy_alph_p_1) {

        int new_cpy_pos = cpy_pos+1;

        // STATE CHANGES: add folds and copy
        fold_r[new_cpy_pos] = 2; //New monomers are always bound to the template. In practice, it can be implemented as an unflippable coil.
        copy_r[new_cpy_pos] = newstate_flag-1; // New monomer in the copy

        // The last added monomer is now a 'proper' coil state (rather than bound to template), insert weights to allow flipping of this coil state
        if (cpy_pos > 3) {
            weight_thresholds[copy_alph_p_1+new_cpy_pos-1] = param_Kfold*exp( GetContigHelix( fold_r[new_cpy_pos-4], fold_r[new_cpy_pos-3], fold_r[new_cpy_pos-2], 2, 0 , 0  )*param_Gu );
            partition_func = partition_func + weight_thresholds[copy_alph_p_1+new_cpy_pos-1]; 
        }
        // update monomer removal weights
        diff = backward_rates_c[ copy_r[new_cpy_pos-1] ][ template_r[new_cpy_pos-1] ]
            [ copy_r[new_cpy_pos] ][ template_r[ new_cpy_pos ] ][ fold_r[new_cpy_pos-1]] - weight_thresholds[0];
        
        weight_thresholds[0] = weight_thresholds[0] + diff;
        partition_func = partition_func + diff;
        
        for ( int i = 0; i < copy_alph_p_1 - 1; i++ ) {
            // update monomer addition weights
            diff = forward_rates_c[ copy_r[new_cpy_pos] ][ template_r[new_cpy_pos] ]
                [ i ][ template_r[new_cpy_pos + 1] ] - weight_thresholds[i+1];

            weight_thresholds[i+1] = weight_thresholds[i+1] + diff;
            partition_func = partition_func + diff;
        };

        return(new_cpy_pos);
    } else {
        int flip_pos = newstate_flag - copy_alph_p_1;
        
        fold_r[flip_pos] =  3 - fold_r[flip_pos];
        
        if ((flip_pos) == cpy_pos - 1) {
            // This is the only possibility that allows for monomer removal rates to change
            // update monomer removal weights
            diff = backward_rates_c[ copy_r[cpy_pos-1] ][ template_r[cpy_pos-1] ]
                [ copy_r[cpy_pos] ][ template_r[ cpy_pos ] ][ fold_r[cpy_pos-1]] - weight_thresholds[0];
        
            weight_thresholds[0] = weight_thresholds[0] + diff;
            partition_func = partition_func + diff;
        }

         if ((flip_pos) == cpy_pos) {
            // error state for debugging only.
            cout << "Flipping at " << cpy_pos << endl;
            throw runtime_error("Flipping of the tip monomer, not allowed as this monomer is bound to the template");
        }

        // Updating the fold
        for( int i = flip_pos-2; i < flip_pos+3; i++ ) {
            if (i == 3) {
                partition_func = partition_func - weight_thresholds[copy_alph_p_1 + i];
                weight_thresholds[copy_alph_p_1 + i] = 0; // start position cannot flip
            } else if (i == cpy_pos) {
                partition_func = partition_func - weight_thresholds[copy_alph_p_1 + i];
                weight_thresholds[copy_alph_p_1 + i] = 0; // The copy position cannot flip
            } else if (fold_r[i] == 1){
                diff = param_Kfold*exp( -GetContigHelix( fold_r[i-3], fold_r[i-2], fold_r[i-1], fold_r[i+1], fold_r[i+2], fold_r[i+3] )
                    *(param_Gh-param_Gu) - param_Gfold  ) - weight_thresholds[copy_alph_p_1 + i]; 
                
                weight_thresholds[copy_alph_p_1 + i] =   weight_thresholds[copy_alph_p_1 + i] + diff; 
                partition_func = partition_func + diff;
            } else if (fold_r[i] == 2) {
                diff = param_Kfold*exp( GetContigHelix( fold_r[i-3], fold_r[i-2], fold_r[i-1], fold_r[i+1], fold_r[i+2], fold_r[i+3] )*param_Gu )
                    - weight_thresholds[copy_alph_p_1 + i];

                weight_thresholds[copy_alph_p_1 + i] =   weight_thresholds[copy_alph_p_1 + i] + diff; 
                partition_func = partition_func + diff;
            } else {
                partition_func = partition_func - weight_thresholds[copy_alph_p_1 + i];
                weight_thresholds[copy_alph_p_1 + i] = 0;
            }
        }

        return(cpy_pos);
    }


}

inline int FoldDrivenSim::UpdateHbondDebug(int newstate_flag, int cpy_pos, int copy_alph_p_1, boostmat<double,5> &backward_rates_c, boostmat<double,4> &forward_rates_c,
            double param_Kfold, double param_Gfold, double param_Gu, double param_Gh, vector<int> &fold_r, vector<int> &copy_r, vector<int> &template_r, 
            vector<double> &weight_thresholds, double &partition_func, int &hydrogen_count){
    
    double diff;

    if (newstate_flag == 0) {
        // A monomer is removed.

        // STATE CHANGES: remove folds and copy 
        fold_r[cpy_pos] = 0;
        copy_r[cpy_pos] = 0;

        // The monomer behind the removed monomer can no longer flip as it is now bound to the copy, set the weight to zero.
        // Update the normalization factor first
        partition_func = partition_func - weight_thresholds[copy_alph_p_1+cpy_pos-1]; //each time a weight is updated, 
            // the normalization factor is also updated. The idea is to avoid having to sum up the whole vector each time
        weight_thresholds[copy_alph_p_1+cpy_pos-1] = 0;

        int new_cpy_pos = cpy_pos-1;
       
        // update monomer removal weights
        diff = backward_rates_c[ copy_r[new_cpy_pos-1] ][ template_r[new_cpy_pos-1] ]
            [ copy_r[new_cpy_pos] ][ template_r[ new_cpy_pos ] ][ fold_r[new_cpy_pos-1]] - weight_thresholds[0];
        
        weight_thresholds[0] = weight_thresholds[0] + diff;
        partition_func = partition_func + diff;

        for ( int i = 0; i < copy_alph_p_1 - 1; i++ ) {
            // update monomer addition weights
            diff = forward_rates_c[ copy_r[new_cpy_pos] ][ template_r[new_cpy_pos] ]
                [ i ][ template_r[new_cpy_pos + 1] ] - weight_thresholds[i+1];

            weight_thresholds[i+1] = weight_thresholds[i+1] + diff;
            partition_func = partition_func + diff;
        };
        
        return(new_cpy_pos);
    } else if (newstate_flag < copy_alph_p_1) {

        int new_cpy_pos = cpy_pos+1;

        // STATE CHANGES: add folds and copy
        fold_r[new_cpy_pos] = 2; //New monomers are always bound to the template. In practice, it can be implemented as an unflippable coil.
        copy_r[new_cpy_pos] = newstate_flag-1; // New monomer in the copy

        // The last added monomer is now a 'proper' coil state (rather than bound to template), insert weights to allow flipping of this coil state
        if (cpy_pos > 3) {
            weight_thresholds[copy_alph_p_1+new_cpy_pos-1] = param_Kfold*exp( GetContigHelix( fold_r[new_cpy_pos-4], fold_r[new_cpy_pos-3], fold_r[new_cpy_pos-2], 2, 0, 0  )*param_Gu );
            partition_func = partition_func + weight_thresholds[copy_alph_p_1+new_cpy_pos-1]; 
        }
        // update monomer removal weights
        diff = backward_rates_c[ copy_r[new_cpy_pos-1] ][ template_r[new_cpy_pos-1] ]
            [ copy_r[new_cpy_pos] ][ template_r[ new_cpy_pos ] ][ fold_r[new_cpy_pos-1]] - weight_thresholds[0];
        
        weight_thresholds[0] = weight_thresholds[0] + diff;
        partition_func = partition_func + diff;
        
        for ( int i = 0; i < copy_alph_p_1 - 1; i++ ) {
            // update monomer addition weights
            diff = forward_rates_c[ copy_r[new_cpy_pos] ][ template_r[new_cpy_pos] ]
                [ i ][ template_r[new_cpy_pos + 1] ] - weight_thresholds[i+1];

            weight_thresholds[i+1] = weight_thresholds[i+1] + diff;
            partition_func = partition_func + diff;
        };

        return(new_cpy_pos);
    } else {
        int flip_pos = newstate_flag - copy_alph_p_1;
        
        // Debugging
        if (fold_r[flip_pos] == 2) {
            cout << "Adding " << GetContigHelix( fold_r[flip_pos-3], fold_r[flip_pos-2], fold_r[flip_pos-1], fold_r[flip_pos+1], fold_r[flip_pos+2], fold_r[flip_pos+3] ) << " to hydrogen" << endl;
            hydrogen_count = hydrogen_count + GetContigHelix( fold_r[flip_pos-3], fold_r[flip_pos-2], fold_r[flip_pos-1], fold_r[flip_pos+1], fold_r[flip_pos+2], fold_r[flip_pos+3] );
        } else {
            cout << "Subtracting " << GetContigHelix( fold_r[flip_pos-3], fold_r[flip_pos-2], fold_r[flip_pos-1], fold_r[flip_pos+1], fold_r[flip_pos+2], fold_r[flip_pos+3] ) << " from hydrogen" << endl;
            hydrogen_count = hydrogen_count - GetContigHelix( fold_r[flip_pos-3], fold_r[flip_pos-2], fold_r[flip_pos-1], fold_r[flip_pos+1], fold_r[flip_pos+2], fold_r[flip_pos+3]);
        }

        fold_r[flip_pos] =  3 - fold_r[flip_pos];
        
        if ((flip_pos) == cpy_pos - 1) {
            // This is the only possibility that allows for monomer removal rates to change
            // update monomer removal weights
            diff = backward_rates_c[ copy_r[cpy_pos-1] ][ template_r[cpy_pos-1] ]
                [ copy_r[cpy_pos] ][ template_r[ cpy_pos ] ][ fold_r[cpy_pos-1]] - weight_thresholds[0];
        
            weight_thresholds[0] = weight_thresholds[0] + diff;
            partition_func = partition_func + diff;
        }

         if ((flip_pos) == cpy_pos) {
            // error state for debugging only.
            throw runtime_error("Flipping of the tip monomer, not allowed as this monomer is bound to the template");
        }

        // Updating the fold
        for( int i = flip_pos-2; i < flip_pos+3; i++ ) {
            if (i == 3) {
                partition_func = partition_func - weight_thresholds[copy_alph_p_1 + i];
                weight_thresholds[copy_alph_p_1 + i] = 0; // start position cannot flip
            } else if (i == cpy_pos) {
                partition_func = partition_func - weight_thresholds[copy_alph_p_1 + i];
                weight_thresholds[copy_alph_p_1 + i] = 0; // The copy position cannot flip
            } else if (fold_r[i] == 1){
                diff = param_Kfold*exp( -GetContigHelix( fold_r[i-3], fold_r[i-2], fold_r[i-1], fold_r[i+1], fold_r[i+2], fold_r[i+3]  )
                    *(param_Gh-param_Gu) - param_Gfold  ) - weight_thresholds[copy_alph_p_1 + i]; 
                
                weight_thresholds[copy_alph_p_1 + i] =   weight_thresholds[copy_alph_p_1 + i] + diff; 
                partition_func = partition_func + diff;
            } else if (fold_r[i] == 2) {
                diff = param_Kfold*exp( GetContigHelix( fold_r[i-3], fold_r[i-2], fold_r[i-1], fold_r[i+1], fold_r[i+2], fold_r[i+3]  )*param_Gu )
                    - weight_thresholds[copy_alph_p_1 + i];

                weight_thresholds[copy_alph_p_1 + i] =   weight_thresholds[copy_alph_p_1 + i] + diff; 
                partition_func = partition_func + diff;
            } else {
                partition_func = partition_func - weight_thresholds[copy_alph_p_1 + i];
                weight_thresholds[copy_alph_p_1 + i] = 0;
            }
        }

        return(cpy_pos);
    }


}

//**UTILITIES**//

inline int FoldDrivenSim::EvalHydrobonds( vector<int> fold ){
    int count = 0;
    for (int i = 3; i < fold.size()-4; i++) {
        if (fold[i + 1] == 1 && fold[i + 2] == 1 && fold[i + 3] == 1) {
            count++;
        }
    }
    return(count);
}

inline int FoldDrivenSim::EvalCoil( vector<int> fold ){
    int count = 0;
    for (int i = 3; i < fold.size()-4; i++) {
        if (fold[i] == 2) {
            count++;
        }
    }
    return(count);
}

double FoldDrivenSim::GetAccuracy() {
    double accuracy = 0;
    for (int i = 4; i < complete_copy.size()-3; i++) {
        accuracy = accuracy + double(template_s[i] == complete_copy[i])/double(complete_copy.size()-7); 
    }
    return(accuracy);
}

void FoldDrivenSim::PrintKineticTensors(){
    cout << "Forward Tensor:" << endl;
    for (int i = 0; i != forward_rates.shape()[0]; i++){
        for (int j = 0; j != forward_rates.shape()[1]; j++){
	        cout << i << j << " " << flush;
            for (int k = 0; k != forward_rates.shape()[2]; k++){
                for (int l = 0; l != forward_rates.shape()[3]; l++){
                    cout << forward_rates [i][j][k][l] << " " << flush;
                }
            }
        cout << endl;
        }       
    }

    cout << "Backward Tensor:" << endl;
    for (int i = 0; i != backward_rates.shape()[0]; i++){
        for (int j = 0; j != backward_rates.shape()[1]; j++){
            cout << i << j << " " << flush;
            for (int k = 0; k != backward_rates.shape()[2]; k++){
                for (int l = 0; l != backward_rates.shape()[3]; l++){
                    cout << backward_rates [i][j][k][l][2] << " " << flush;
                }
            }
            cout << endl;
        }
    }

}

#endif
