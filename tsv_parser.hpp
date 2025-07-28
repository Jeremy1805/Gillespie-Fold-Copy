#ifndef _TsvParser_H
#define _TsvParser_H

#include "boost/multi_array.hpp"

template<class T, int N>
using boostmat = typename boost::multi_array<T,N>;

typedef boost::multi_array<int, 2> array_type;
typedef boost::multi_array_types::index_range range_t;

using namespace std;

void PrintMatrix(boostmat<double,2> mat_print);

template<class T>
void PrintVector(vector<T> vec_print);

// Function for getting back the boostmatrix from the string

template<class T>
boostmat<T,2> Str2Mat(int row, int col, string mat_string);

// Template declaration

template<class T>
vector<T> Str2Vec(string vec_string);

// Function for getting an int vector from a string

template<>
vector<int> Str2Vec<int>(string vec_string);

// Get an unsigned int vector from a semicolon-delimited string

template<>
vector<unsigned int> Str2Vec<unsigned int>(string vec_string); 

// Get a vector of doubles from a semicolon-delimited string

template<>
vector<double> Str2Vec<double>(string vec_string);

// Function for getting unordered map from a string 
// that defines an unordered map

class TsvParser {
    
    public:

    /*0*/ string batch_id;            // One ID for the whole csv file
    /*1*/ string data_id;             // One ID for each data
    /*2*/ string template_label;  // Fast processing of Templates
    /*3*/ int copy_alph_size;     // Size of the copy alphabet
    /*4*/ int template_alph_size; // Size of the template alphabet
    /*5*/ int template_length;
    
    //FoldDrivenSim Class Inputs
    /*6*/ vector<int> initial_copy;
    /*7*/ vector<int> initial_template;
    /*8*/ vector<int> initial_fold;
    /*9*/ double G_pol;
    /*10*/ double k_copy;
    /*11*/ double k_fold;
    /*12*/ double G_fold;
    /*13*/ double G_u;
    /*14*/ double G_h;
    /*15*/ boostmat<double,2> energy_matrix;
    /*16*/ int kinetic_flag;
    /*17*/ boostmat<double,2> fold_map;

    /*18*/ vector<unsigned int> random_seed;
    /*19*/ int simulation_repeat;
    /*20*/ string output_file_name;
    /*21*/ bool debug;
    
    // Output only fieldsFoldDrivenSim
    /*22*/ vector<vector<int>> output_copy_list;
    /*23*/ vector<vector<int>> output_folds_list;
    /*24*/ vector<double> accuracy_list;
    /*25*/ double mean_accuracy;
    
    TsvParser();
    TsvParser(const TsvParser &tsv_parser);
    TsvParser& operator=(const TsvParser &tsv_parser);
    TsvParser( string ); // Constructor for reading from a file 
    TsvParser( string batch_id, string data_id, string template_in, vector<int> initial_template_in, double G_pol_in, double k_copy_in,
        double k_fold_in, double g_fold_in, double g_u_in, double g_h_in, boostmat<double,2> energy_matrix_in, int kinetic_flag_in, boostmat<double,2> fold_map_in,
        int simulation_repeat_in, string output_file_in); // Constructor for creating a parameter set that generates its own random numbers, and also the single-base initial copy
    string Stringify(); // Generates a string for output
    string StringifyOutput();
    string StringifySummary();
    string Header();
    string OutputHeader();
    string SummaryHeader();
    void Print();
};

#endif
