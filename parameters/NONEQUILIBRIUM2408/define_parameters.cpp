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
#include <eigen3/Eigen/Eigenvalues> // header file
#include "tsv_parser.cpp"

template<class T, int N>
using boostmat = typename boost::multi_array<T,N>;

typedef boost::multi_array<int, 2> array_type;
typedef boost::multi_array_types::index_range range_t;

using namespace std;

int main(){
    /* Generate_batch_id */
    string batch = "NONEQUILIBRIUM2408";
    // For a G_pol calculated at the equilibrium energy -ln(lambda), 
       // how does efficicency change with increasing speed of folding k_fold
       // for various fixed levels of G_h and G_fold 

    vector<TsvParser> tsv_parser_vec;
    
    vector<int> template_s(507,1);
    for (int i = 0; i < 3; i++) {
        template_s[i] = 0;
        template_s[template_s.size()-1-i] = 0;
    }
    
    vector<int> copy(507,0);
    copy[3] = 1; 
    vector<int> fold(507,0);
    fold[3] = 1; 
    
    int repeats = 50;
    
    // Fixing G_k_tt at 2
    boostmat<double,2> energy_matrix;
    energy_matrix.resize(boost::extents[3][3]);
    energy_matrix[1][1] = 4;
    energy_matrix[2][2] = 4;
    
    string output_file = batch + ".results";
    int id_counter = 0;

    for (double G_fold = 0.2; G_fold < 1.1; G_fold = G_fold + 0.2 ){
        for (double G_h = 0.2; G_h < 1.1; G_h = G_h + 0.2){
            double w = exp(G_h+G_fold);
            double v = exp(G_fold);
            double u = 1;

            Eigen::Matrix<double, 3, 3>  lifson_roig_transfer; // declare a real (double) 2x2 matrix
            lifson_roig_transfer << w, v, 0, 0, 0, u, v, v, u; // defined the matrix A

            Eigen::EigenSolver<Eigen::Matrix<double, 3,3> > eigen_lifson(lifson_roig_transfer); // the instance s(A) includes the eigensystem
            
            for ( double G_pol_base = -0.6; G_pol_base < 0.05; G_pol_base = G_pol_base + 0.1) {
                double G_pol = G_pol_base - log(eigen_lifson.eigenvalues().real().maxCoeff());
                
                for (int k_fold = 2; k_fold < 21; k_fold = k_fold+2) { 
                    TsvParser tsv_parser( batch, to_string(id_counter), "111", template_s, G_pol, 1,
                        k_fold, G_fold, 0, G_h, energy_matrix, 0, energy_matrix, repeats, output_file);
                    tsv_parser_vec.push_back(tsv_parser);
                    id_counter++;
                }

                TsvParser tsv_parser( batch, to_string(id_counter), "111", template_s, G_pol, 1,
                        100, G_fold, 0, G_h, energy_matrix, 0, energy_matrix, repeats, output_file);
                    tsv_parser_vec.push_back(tsv_parser);
                id_counter++;
            }
            
            for ( double G_pol_base = 0.2; G_pol_base < 1; G_pol_base = G_pol_base + 0.2) {
                double G_pol = G_pol_base - log(eigen_lifson.eigenvalues().real().maxCoeff());
                
                for (int k_fold = 2; k_fold < 21; k_fold = k_fold+2) { 
                    TsvParser tsv_parser( batch, to_string(id_counter), "111", template_s, G_pol, 1,
                        k_fold, G_fold, 0, G_h, energy_matrix, 0, energy_matrix, repeats, output_file);
                    tsv_parser_vec.push_back(tsv_parser);
                    id_counter++;
                }

                TsvParser tsv_parser( batch, to_string(id_counter), "111", template_s, G_pol, 1,
                        100, G_fold, 0, G_h, energy_matrix, 0, energy_matrix, repeats, output_file);
                    tsv_parser_vec.push_back(tsv_parser);
                id_counter++;
            }

        }
    }
    
    for ( double G_pol_base = -0.6; G_pol_base < 1; G_pol_base = G_pol_base + 0.1) {
        double G_pol = G_pol_base;        
        for (int k_fold = 2; k_fold < 21; k_fold = k_fold+2) { 
            TsvParser tsv_parser( batch, to_string(id_counter), "111", template_s, G_pol, 1,
                    0, 0, 0, 0, energy_matrix, 0, energy_matrix, repeats, output_file);
            tsv_parser_vec.push_back(tsv_parser);
            id_counter++;
        }
    }
    
    std::ofstream tsv_output;
    tsv_output.open(batch + ".input.tsv");
    tsv_output << TsvParser().Header() << endl;
    for ( TsvParser single_data : tsv_parser_vec){
        tsv_output << single_data.Stringify() << endl;
    }
}
