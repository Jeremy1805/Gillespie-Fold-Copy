#ifndef _TsvParse
#define _TsvParse

#include <iostream>
#include <fstream>
#include <algorithm>
#include <tuple>
#include <string>
#include <ctype.h>
#include <random>
#include <iomanip>
#include <unordered_map>

#include "boost/multi_array.hpp"
#include <boost/algorithm/string.hpp>
#include "boost/lexical_cast.hpp"

#include "tsv_parser.hpp"

template<class T, int N>
using boostmat = typename boost::multi_array<T,N>;

typedef boost::multi_array<int, 2> array_type;
typedef boost::multi_array_types::index_range range_t;

using namespace std;

// Function for printing a matrix

void PrintMatrix(boostmat<double,2> mat_print){
    for(int x = 0; x < mat_print.shape()[0]; x++) {
        for(int y = 0; y < mat_print.shape()[1]; y++) {
            cout << mat_print[x][y] << " " << flush;
        }
	cout << endl;
    }
}

template<class K, class V>
void PrintUnorderedMap(unordered_map<K,V> map_print){
    for (pair<K, V> element : map_print)
    {
        cout << element.first << ": " << element.second << endl;
    }
}

template<class T>
void PrintVector(vector<T> vec_print){
    for (T element : vec_print)
    {
        cout << to_string(element) << " " << flush;
    }
    cout << endl;
}

// Functor to check that only allowed values 
// are present in the energy matrix

struct CheckMatrixString
{
  bool stringIsCorrect = true;

  void operator() (char ch)
  {
    if(stringIsCorrect && !(((ch >= '0' )&&(ch <= '9')) || 
            (ch == ',') || (ch == ';') || (ch == '.') || (ch == '-'))){
	    cout << "UNEXPECTED CHARACTER: " << ch << endl;
            stringIsCorrect = false;
    }
  }
};

// Function for getting back the boostmatrix from the string

template<class T>
boostmat<T,2> Str2Mat(int row, int col, string mat_string){
    boostmat<T,2> out_matrix{boost::extents[row][col]};
    
    mat_string.erase(remove_if(mat_string.begin(), 
    	mat_string.end(), ::isspace), mat_string.end());
    
    CheckMatrixString checker;
    for_each(mat_string.begin(), mat_string.end(), ref(checker));
    
    if(!checker.stringIsCorrect) {
        cout << "Unexpected character in function" <<
                " Str2Mat! Ensure that the Matrix only" <<
	        " has characters '0'-'9','.',',','-' and ';'!"  << endl;
        cout << "Offending string is " << mat_string << endl;
        throw;
    }

    // Delimit the string to get vector of string elements
    vector<string> string_elements;
    //remove trailing delimiters first
    boost::trim_if(mat_string, boost::is_any_of(",;")); 
    boost::split(string_elements, mat_string, boost::is_any_of(",;"));
    
    // Convert from str to double
    vector<double> double_elements(string_elements.size());
    transform(string_elements.begin(), string_elements.end(), 
	    double_elements.begin(),
            [](string const& val) {return stod(val);});
    
    // Throw error if double_elements has more elements than expected
    if (double_elements.size() != row*col) {
        cout << "Size mismatch between expected" <<
	        " and actual matrix size in Str2Mat" << endl;
	throw;
    }
    
    // Assign elements to the boost array
    
    int counter = 0;
    
    for(int x = 0; x < row; x++) {
        for(int y = 0; y < col; y++) {
            out_matrix[x][y] = double_elements[counter];
	    counter++;
        }
    }
    
    return out_matrix;
}

// Template declaration

template<class T>
vector<T> Str2Vec(string vec_string);

// Function for getting an int vector from a string

template<>
vector<int> Str2Vec<int>(string vec_string) {
    vector<int> out_vec;
    // Remove whitespace	
    vec_string.erase(remove_if(vec_string.begin(), 
        vec_string.end(), ::isspace), vec_string.end());

    // Initial delimiting to get vector of strings
    vector<string> vec_string_split;
    
    // remove trailing delimiters first
    boost::trim_if(vec_string, boost::is_any_of(";")); 
    boost::split(vec_string_split, vec_string, boost::is_any_of(";"));
    
    // Within each key-value pair, separate keys from values
    for ( string element: vec_string_split ) {
	out_vec.push_back(stoi(element));
    }
    
    return out_vec;
}

// Get an unsigned int vector from a semicolon-delimited string

template<>
vector<unsigned int> Str2Vec<unsigned int>(string vec_string) {
    vector<unsigned int> out_vec;
    // Remove whitespace	
    vec_string.erase(remove_if(vec_string.begin(), 
        vec_string.end(), ::isspace), vec_string.end());

    // Initial delimiting to get vector of strings
    vector<string> vec_string_split;
    
    // remove trailing delimiters first
    boost::trim_if(vec_string, boost::is_any_of(";")); 
    boost::split(vec_string_split, vec_string, boost::is_any_of(";"));
    
    // Within each key-value pair, separate keys from values
    for ( string element: vec_string_split ) {
	out_vec.push_back(stoul(element));
    }
    
    return out_vec;
}

// Get a vector of doubles from a semicolon-delimited string

template<>
vector<double> Str2Vec<double>(string vec_string) {
    vector<double> out_vec;
    // Remove whitespace	
    vec_string.erase(remove_if(vec_string.begin(), 
        vec_string.end(), ::isspace), vec_string.end());

    // Initial delimiting to get vector of strings
    vector<string> vec_string_split;
    
    // remove trailing delimiters first
    boost::trim_if(vec_string, boost::is_any_of(";")); 
    boost::split(vec_string_split, vec_string, boost::is_any_of(";"));
    
    // Within each key-value pair, separate keys from values
    for ( string element: vec_string_split ) {
	out_vec.push_back(stod(element));
    }
    
    return out_vec;
}

// Function for getting unordered map from a string 
// that defines an unordered map
template<class T>
string Mat2Str(boostmat<T,2> mat2sav){
    string matrix_string = "";
    for(int x = 0; x < mat2sav.shape()[0]; x++) {
        for(int y = 0; y < mat2sav.shape()[1]; y++) {
            matrix_string = matrix_string + to_string(mat2sav[x][y]);
	    if (y == mat2sav.shape()[1]-1){
		if (x != mat2sav.shape()[0]-1){
		    // Distinction between ',' and ';' is only
		    // for human readability. The program treats
		    // them the same
                    matrix_string = matrix_string + ";";
		}
	    } else {
		matrix_string = matrix_string + ",";
	    }
        }
    }
    return matrix_string;
}

string get_batch_id(){
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    
    std::ostringstream oss;
    oss << std::put_time(&tm, "%d%m%Y%H%M%S");
    auto str = oss.str();
    
    return str;
}

// Function for converting a vector
// to a semi-colon delimited string
template<class K>
string Vec2Str(vector<K> data_vec){
    string vec_string = "";
    
    // Iterate over each key-value pair and turn it to a string
    for (K data_element : data_vec)
    {
        vec_string = vec_string + to_string(data_element);
        vec_string = vec_string + ";"; 
    }
    
    // Remove trailing ';'
    boost::trim_if(vec_string, boost::is_any_of(";"));

    return vec_string;
}

// Create an instance of the tsv parser class - reads a single row of the tsv and stores that row
// as an instance of the class
TsvParser::TsvParser() {}
TsvParser::TsvParser(string input_string) {
    vector<string> tsv_cells;
    split(tsv_cells, input_string, boost::is_any_of("\t"));
    // Check that the number of columns is the expected number of columns
    if (tsv_cells.size() != 22){
        cout << "ERROR: INCORRECT COLUMN SIZE OF TSV FILE!" << endl;
        throw;    
    }
    
    batch_id = tsv_cells[0];            // One ID for the whole csv file
    data_id = tsv_cells[1];             // One ID for each data
    template_label = tsv_cells[2];  // Fast processing of Templates
    
    copy_alph_size = stoi(tsv_cells[3]);     // Size of the copy alphabet
    template_alph_size = stoi(tsv_cells[4]); // Size of the template alphabet
    template_length = stoi(tsv_cells[5]);
    
    initial_copy = Str2Vec<int>(tsv_cells[6]);
    initial_template = Str2Vec<int>(tsv_cells[7]);
    initial_fold = Str2Vec<int>(tsv_cells[8]);

    G_pol = stod(tsv_cells[9]);
    k_copy = stod(tsv_cells[10]);
    k_fold = stod(tsv_cells[11]);
    G_fold = stod(tsv_cells[12]);
    G_u = stod(tsv_cells[13]);
    G_h = stod(tsv_cells[14]);

    energy_matrix.resize(boost::extents[copy_alph_size][template_alph_size]);
    energy_matrix = Str2Mat<double>(copy_alph_size, template_alph_size, tsv_cells[15]);
    
    kinetic_flag = stoi(tsv_cells[16]);

    fold_map.resize(boost::extents[copy_alph_size][3]);
    fold_map = Str2Mat<double>(copy_alph_size,3,tsv_cells[17]);

    random_seed = Str2Vec<unsigned int>(tsv_cells[18]);;
    simulation_repeat = stoi(tsv_cells[19]);
    output_file_name = tsv_cells[20];
    debug = boost::lexical_cast<bool>(tsv_cells[21]);
}

TsvParser::TsvParser(const TsvParser &tsv_parser) {
    batch_id = tsv_parser.batch_id;            // One ID for the whole csv file
    data_id = tsv_parser.data_id;             // One ID for each data 
    template_label = tsv_parser.template_label;      // Label describing how the template was created
    
    // Sizes and lengths
    copy_alph_size = tsv_parser.copy_alph_size;     // Size of the copy alphabet
    template_alph_size = tsv_parser.template_alph_size; // Size of the template alphabet
    template_length = tsv_parser.template_length;
    
    initial_copy = tsv_parser.initial_copy;
    initial_template = tsv_parser.initial_template;
    initial_fold = tsv_parser.initial_fold;

    // Short parameters
    G_pol = tsv_parser.G_pol;
    k_copy = tsv_parser.k_copy;
    k_fold = tsv_parser.k_fold;
    G_fold = tsv_parser.G_fold;
    G_u = tsv_parser.G_u;
    G_h = tsv_parser.G_h;
    
    energy_matrix.resize(boost::extents[copy_alph_size][template_alph_size]);
    energy_matrix = tsv_parser.energy_matrix;
    
    kinetic_flag = tsv_parser.kinetic_flag;

    fold_map.resize(boost::extents[copy_alph_size][3]);
    fold_map = tsv_parser.fold_map;
    
    // Miscellaneous
    random_seed = tsv_parser.random_seed;
    simulation_repeat = tsv_parser.simulation_repeat;
    output_file_name = tsv_parser.output_file_name;
    debug = tsv_parser.debug;

    //Outputs
    output_copy_list = tsv_parser.output_copy_list;
    output_folds_list = tsv_parser.output_folds_list;
    accuracy_list = tsv_parser.accuracy_list;
    mean_accuracy = tsv_parser.mean_accuracy;
}

TsvParser& TsvParser::operator=(const TsvParser &tsv_parser)
{
    batch_id = tsv_parser.batch_id;            // One ID for the whole csv file
    data_id = tsv_parser.data_id;             // One ID for each data 
    template_label = tsv_parser.template_label;      // Label describing how the template was created
    
    // Sizes and lengths
    copy_alph_size = tsv_parser.copy_alph_size;     // Size of the copy alphabet
    template_alph_size = tsv_parser.template_alph_size; // Size of the template alphabet
    template_length = tsv_parser.template_length;
    
    initial_copy = tsv_parser.initial_copy;
    initial_template = tsv_parser.initial_template;
    initial_fold = tsv_parser.initial_fold;

    // Short parameters
    G_pol = tsv_parser.G_pol;
    k_copy = tsv_parser.k_copy;
    k_fold = tsv_parser.k_fold;
    G_fold = tsv_parser.G_fold;
    G_u = tsv_parser.G_u;
    G_h = tsv_parser.G_h;
    
    energy_matrix.resize(boost::extents[copy_alph_size][template_alph_size]);
    energy_matrix = tsv_parser.energy_matrix;
    
    kinetic_flag = tsv_parser.kinetic_flag;

    fold_map.resize(boost::extents[copy_alph_size][3]);
    fold_map = tsv_parser.fold_map;
    
    // Miscellaneous
    random_seed = tsv_parser.random_seed;
    simulation_repeat = tsv_parser.simulation_repeat;
    output_file_name = tsv_parser.output_file_name;
    debug = tsv_parser.debug;

    //Outputs
    output_copy_list = tsv_parser.output_copy_list;
    output_folds_list = tsv_parser.output_folds_list;
    accuracy_list = tsv_parser.accuracy_list;
    mean_accuracy = tsv_parser.mean_accuracy;

    // return the existing object so we can chain this operator
    return *this;
}

TsvParser::TsvParser( string batch_id_in, string data_id_in, string template_label_in, vector<int> initial_template_in, double G_pol_in, double k_copy_in,
        double k_fold_in, double g_fold_in, double g_u_in, double g_h_in, boostmat<double,2> energy_matrix_in, int kinetic_flag_in, boostmat<double,2> fold_map_in,
        int simulation_repeat_in, string output_file_in){
    
    /*0*/ batch_id = batch_id_in;            // One ID for the whole csv file
    /*1*/ data_id = data_id_in;             // One ID for each data
    /*2*/ template_label = template_label_in;  // Fast processing of Templates
    /*3*/ copy_alph_size = energy_matrix_in.shape()[0];     // Size of the copy alphabet
    /*4*/ template_alph_size = energy_matrix_in.shape()[1]; // Size of the template alphabet
    /*5*/ template_length = initial_template_in.size();
    
    //FoldDrivenSim Class Inputs
    /*6*/ initial_copy = vector<int>(template_length,0);
                initial_copy[3] = initial_template_in[3]; 
    /*7*/ initial_template = initial_template_in;
    /*8*/ initial_fold = vector<int>(template_length,0);
                initial_fold[3] = 2;

    /*9*/  G_pol = G_pol_in;
    /*10*/ k_copy = k_copy_in;
    /*11*/ k_fold = k_fold_in;
    /*12*/ G_fold = g_fold_in;
    /*13*/ G_u = g_u_in;
    /*14*/ G_h = g_h_in;

    /*15*/  energy_matrix.resize(boost::extents[copy_alph_size][template_alph_size]);
                energy_matrix = energy_matrix_in;
    
    /*16*/  kinetic_flag = 0;
    
    /*17*/  fold_map.resize(boost::extents[copy_alph_size][3]);
                fold_map = fold_map_in;

    /*18*/  random_device device;
                for (int i = 0; i < simulation_repeat_in; i++) {
                    random_seed.push_back(device());
                }
    /*19*/ simulation_repeat = simulation_repeat_in;
    /*20*/ output_file_name = output_file_in;
    
    /*21*/ debug = 0;

}

string TsvParser::Header(){
    return "BatchID\tParamSetID\tTemplateLabel\t"
    "CopyAlphabetSize\tTemplateAlphabetSize"
    "\tTemplateLength\tInitCopy\tTemplate\tInitFold\tG_pol\t"
    "k_copy\tk_fold\tG_fold\tG_u\tG_h\tEnergyMatrix\t"
    "KineticFlag\tFoldMap" 
    "\tRandomSeeds\tSimulationRepeats\tOutputFile\tDebug";
}

string TsvParser::OutputHeader(){
    return "BatchID\tParamSetID\tTemplateLabel\t"
    "CopyAlphabetSize\tTemplateAlphabetSize"
    "\tTemplateLength\tInitCopy\tTemplate\tInitFold\tG_pol\t"
    "k_copy\tk_fold\tG_fold\tG_u\tG_h\tEnergyMatrix\t"
    "KineticFlag\tFoldMap" 
    "\tRandomSeeds\tSimulationRepeats\tOutputFile\tDebug\tCopyFinal\tFoldFinal\tAccuracy";
}

string TsvParser::SummaryHeader(){
    return "BatchID\tParamSetID\tTemplateLabel\t"
    "CopyAlphabetSize\tTemplateAlphabetSize"
    "\tTemplateLength\tInitCopy\tTemplate\tInitFold\tG_pol\t"
    "k_copy\tk_fold\tG_fold\tG_u\tG_h\tEnergyMatrix\t"
    "KineticFlag\tFoldMap" 
    "\tRandomSeeds\tSimulationRepeats\tOutputFile\tDebug\tMeanAccuracy";
}

string TsvParser::Stringify(){
    return 
    /*0*/ batch_id + "\t" + 
    /*1*/ data_id + "\t" +
    /*2*/ template_label + "\t" +
    /*3*/ to_string(copy_alph_size) + "\t" + 
    /*4*/ to_string(template_alph_size) + "\t" + 
    /*5*/ to_string(template_length) + "\t" + 
    /*6*/ Vec2Str<int>(initial_copy) + "\t" +
    /*7*/ Vec2Str<int>(initial_template) + "\t" +
    /*8*/ Vec2Str<int>(initial_fold) + "\t" +
    /*9*/ to_string(G_pol) + "\t" +
    /*10*/ to_string(k_copy) + "\t" +
    /*11*/ to_string(k_fold) + "\t" +
    /*12*/ to_string(G_fold) + "\t" +
    /*13*/ to_string(G_u) + "\t" +
    /*14*/ to_string(G_h) + "\t" +
    /*15*/ Mat2Str<double>(energy_matrix) + "\t" + 
    /*16*/ to_string(int(kinetic_flag)) + "\t" +
    /*17*/ Mat2Str(fold_map) + "\t" + 
    /*18*/ Vec2Str<unsigned int>(random_seed) + "\t" +
    /*19*/ to_string(simulation_repeat) + "\t" + 
    /*20*/ output_file_name + "\t" + 
    /*21*/ to_string(debug);
}

string TsvParser::StringifyOutput(){
    if (output_copy_list.size() != random_seed.size()) {
        throw runtime_error("Random seed and output size mismatch!");
    }

    if (output_copy_list.size() != output_folds_list.size() || output_copy_list.size() != accuracy_list.size() )  {
        string error = string("Mismatching copy, fold and accuracy sizes!") + string(" copy: ") 
            + to_string(int(output_copy_list.size())) + string(" folds: ") + to_string(int(output_folds_list.size())) +  string(" accuracy: ") 
            + to_string(int(accuracy_list.size()));
        throw runtime_error(error);
    }
    string repeat_string = 
    /*0*/ batch_id + "\t" + 
    /*1*/ data_id + "\t" +
    /*2*/ template_label + "\t" +
    /*3*/ to_string(copy_alph_size) + "\t" + 
    /*4*/ to_string(template_alph_size) + "\t" + 
    /*5*/ to_string(template_length) + "\t" + 
    /*6*/ Vec2Str<int>(initial_copy) + "\t" +
    /*7*/ Vec2Str<int>(initial_template) + "\t" +
    /*8*/ Vec2Str<int>(initial_fold) + "\t" +
    /*9*/ to_string(G_pol) + "\t" +
    /*10*/ to_string(k_copy) + "\t" +
    /*11*/ to_string(k_fold) + "\t" +
    /*12*/ to_string(G_fold) + "\t" +
    /*13*/ to_string(G_u) + "\t" +
    /*14*/ to_string(G_h) + "\t" +
    /*15*/ Mat2Str<double>(energy_matrix) + "\t" + 
    /*16*/ to_string(int(kinetic_flag)) + "\t" +
    /*17*/ Mat2Str(fold_map);
    string out_string = "";
    for (int i = 0; i < output_copy_list.size(); i++) {
        out_string = out_string + repeat_string + "\t" + 
            to_string(random_seed[i]) + "\t" +
            to_string(simulation_repeat) + "\t" + 
            output_file_name + "\t" + 
            to_string(debug) + "\t" + 
            Vec2Str<int>(output_copy_list[i]) + "\t" + 
            Vec2Str<int>(output_folds_list[i]) + "\t" + 
            to_string(accuracy_list[i]) + "\n";
    }
    return(out_string);
}

string TsvParser::StringifySummary(){
    return 
    /*0*/ batch_id + "\t" + 
    /*1*/ data_id + "\t" +
    /*2*/ template_label + "\t" +
    /*3*/ to_string(copy_alph_size) + "\t" + 
    /*4*/ to_string(template_alph_size) + "\t" + 
    /*5*/ to_string(template_length) + "\t" + 
    /*6*/ Vec2Str<int>(initial_copy) + "\t" +
    /*7*/ Vec2Str<int>(initial_template) + "\t" +
    /*8*/ Vec2Str<int>(initial_fold) + "\t" +
    /*9*/ to_string(G_pol) + "\t" +
    /*10*/ to_string(k_copy) + "\t" +
    /*11*/ to_string(k_fold) + "\t" +
    /*12*/ to_string(G_fold) + "\t" +
    /*13*/ to_string(G_u) + "\t" +
    /*14*/ to_string(G_h) + "\t" +
    /*15*/ Mat2Str<double>(energy_matrix) + "\t" + 
    /*16*/ to_string(int(kinetic_flag)) + "\t" +
    /*17*/ Mat2Str(fold_map) + "\t" + 
    /*18*/ Vec2Str<unsigned int>(random_seed) + "\t" +
    /*19*/ to_string(simulation_repeat) + "\t" + 
    /*20*/ output_file_name + "\t" + 
    /*21*/ to_string(debug) + "\t" + 
    /*22*/ to_string(mean_accuracy);

}

void TsvParser::Print() {
    cout << "BatchID: "<< batch_id << "\t" << flush;
    cout << "DataID: " << data_id << "\t" << flush;
    cout << "TemplateLabel: " << template_label << endl;
    
    cout << "CopyAlphabetSize: "<< to_string(copy_alph_size) << "\t" << flush;
    cout << "TemplateAlphabetSize: " << to_string(template_alph_size) << "\t" << flush;
    cout << "TemplateLength: " << to_string(template_length) << "\t" << flush;

    cout << "CopySequence: " << flush;
    PrintVector(initial_copy);
    cout << "TemplateSequence: " << flush;
    PrintVector(initial_template);
    cout << "FoldSequence: " << flush;
    PrintVector(initial_fold);

    cout << "G_pol: " << to_string(G_pol) << endl; 
    cout << "k_copy: " << to_string(k_copy) << endl; 
    cout << "k_fold: " << to_string(k_fold) << endl; 
    cout << "G_fold: " << to_string(G_fold) << endl; 
    cout << "G_u: " << to_string(G_u) << endl; 
    cout << "G_h: " << to_string(G_h) << endl; 
    cout << "EnergyMatrix: " << endl;
    PrintMatrix(energy_matrix);
    
    cout << "kinetic_flag: " << to_string(int(kinetic_flag)) << endl; 
    cout << "FoldMap: " << endl;
    PrintMatrix(fold_map);

    cout << "RandomSeed: " << flush;
    PrintVector(random_seed);
    cout << "Repeats: "<< to_string(simulation_repeat) << "\t" << flush;
    cout << "OutputFileName: " << output_file_name << "\t" << flush;
    cout << "Debug: " << to_string(debug) << "\t" << flush; 

    cout << endl;
}

void Debug() {
    /*Minimal Labels*/
    /*0*/ string batch_id = get_batch_id();
            
    /*1*/ string data_id = to_string(0);             // One ID for each data
    /*2*/ string template_label = "all_1";  // Fast processing of Templates
    
    //FoldDrivenSim Class Inputs
    /*7*/ vector<int> initial_template(1007,1);
    for (int i = 0; i < 3; i++) {
        initial_template[i] = 0;
        initial_template[ initial_template.size()-1-i ] = 0;
    }

    /*8*/ double G_pol = 0.3;
    /*9*/ double k_copy = 1;
    /*10*/ double k_fold = 10;
    /*11*/ double G_fold = 0.2;
    /*12*/ double G_u = 0;
    /*13*/ double G_h = 0.3;
    /*14*/ boostmat<double,2> energy_matrix;
        energy_matrix.resize(boost::extents[3][3]);
        energy_matrix[1][1] = 1;
        energy_matrix[2][1] = 1;
        energy_matrix[1][2] = 1;
        energy_matrix[2][2] = 1;

    /*15*/ int kinetic_flag = 0;
    /*16*/ boostmat<double,2> fold_map;
        fold_map.resize(boost::extents[3][3]);

    /*18*/ int simulation_repeat = 50;
    /*19*/ string output_file_name = "output.txt";
    
    TsvParser parameter_generator(batch_id, data_id, template_label, initial_template, G_pol, k_copy,
        k_fold, G_fold, G_u, G_h, energy_matrix, kinetic_flag, fold_map,
        simulation_repeat, output_file_name);
    
    std::ofstream tsv_output;
    tsv_output.open(batch_id + ".input.tsv");
    tsv_output << parameter_generator.Header() << endl;
    tsv_output << parameter_generator.Stringify() << endl;
    tsv_output.close();
    
    std::ifstream tsv_input;
    tsv_input.open(batch_id + ".input.tsv");
    
    string first_line;
    string line;
    
    cout << "Checking first line for header.." << endl;

    // Check first line 
    //
    getline(tsv_input, first_line);
    getline(tsv_input,line);
    
    tsv_input.close();

    TsvParser output_tsv(line);
    output_tsv.Print();
}

#endif
