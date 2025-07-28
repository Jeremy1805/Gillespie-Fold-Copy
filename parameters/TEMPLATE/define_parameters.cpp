#include <iostream>
#include <algorithm>
#include <tuple>
#include <string>
#include <ctype.h>
#include <unordered_map>
#include <ctime>
#include <iomanip>
#include <random>
#include <fstream>
#include <cmath>

#include "boost/multi_array.hpp"
#include <boost/algorithm/string.hpp>
#include "tsv_parser.cpp"

template<class T, int N>
using boostmat = typename boost::multi_array<T,N>;

typedef boost::multi_array<int, 2> array_type;
typedef boost::multi_array_types::index_range range_t;

using namespace std;

int main(){
    /* Generate_batch_id */
    string batch = get_batch_id();
    
    vector<TsvParser> tsv_parser_vec;
    
    vector<int> template_s(1007,1);
    for (int i = 0; i < 3; i++) {
        template_s[i] = 0;
        template_s[template_s.size()-1-i] = 0;
    }
    
    vector<int> copy(1007,0);
    copy[3] = 1; 
    vector<int> fold(1007,0);
    fold[3] = 1; 
    
    int repeats = 50;
    
    boostmat<double,2> energy_matrix;
    energy_matrix.resize(boost::extents[3][3]);

    string output_file = batch + ".results";
    int id_counter = 0;
    double lowest_Gpol = -0.2;
    for (double G_k_tt = 0.5; G_k_tt < 8; G_k_tt = G_k_tt + 4.5 ){
        energy_matrix[1][1] = G_k_tt;
        energy_matrix[2][2] = G_k_tt;
        for (double G_pol = lowest_Gpol; G_pol < 0.3; G_pol = G_pol + 0.2){
            TsvParser tsv_no_fold( batch, to_string(id_counter), "111", template_s, G_pol, 1,
                0, 0, 0, 0, energy_matrix, 0, energy_matrix, repeats, output_file);
            id_counter++;
            //g_fold g_h lambda
            // 0.2 1 3.52
            // 0.4 1 4.24
            // 0.6 1 5.12
            // 0.2 0.5 2.56
            // 0.5 0.5 3.21
            // 1   0.5 4.926

            // -true_G_pol = -param_G_pol - ln(lambda) 
            // param G_pol = true_G_pol - ln(lambda)
            TsvParser tsv_fold_1( batch, to_string(id_counter), "111", template_s, G_pol - log(3.52), 1,
                100, 0.2, 0, 1, energy_matrix, 0, energy_matrix, repeats, output_file);
            id_counter++;

            TsvParser tsv_fold_2( batch, to_string(id_counter), "111", template_s, G_pol - log(4.24), 1,
                100, 0.4, 0, 1, energy_matrix, 0, energy_matrix, repeats, output_file);
            id_counter++;

            TsvParser tsv_fold_3( batch, to_string(id_counter), "111", template_s, G_pol - log(5.12), 1,
                100, 0.6, 0, 1, energy_matrix, 0, energy_matrix, repeats, output_file);
            id_counter++;

            TsvParser tsv_fold_4( batch, to_string(id_counter), "111", template_s, G_pol - log(2.56), 1,
                100, 0.2, 0, 0.5, energy_matrix, 0, energy_matrix, repeats, output_file);
            id_counter++;

            TsvParser tsv_fold_5( batch, to_string(id_counter), "111", template_s, G_pol - log(3.21), 1,
                100, 0.5, 0, 0.5, energy_matrix, 0, energy_matrix, repeats, output_file);
            id_counter++;

            TsvParser tsv_fold_6( batch, to_string(id_counter), "111", template_s, G_pol - log(4.93), 1,
                100, 1, 0, 0.5, energy_matrix, 0, energy_matrix, repeats, output_file);
            id_counter++;

            tsv_parser_vec.push_back(tsv_no_fold);
            tsv_parser_vec.push_back(tsv_fold_1);
            tsv_parser_vec.push_back(tsv_fold_2);
            tsv_parser_vec.push_back(tsv_fold_3);
            tsv_parser_vec.push_back(tsv_fold_4);
            tsv_parser_vec.push_back(tsv_fold_5);
            tsv_parser_vec.push_back(tsv_fold_6);
        }
    }

    
    std::ofstream tsv_output;
    tsv_output.open(batch + ".input.tsv");
    tsv_output << TsvParser().Header() << endl;
    for ( TsvParser single_data : tsv_parser_vec){
        tsv_output << single_data.Stringify() << endl;
    }
}
