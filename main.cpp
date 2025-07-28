// Your First C++ Program

#include <iostream>
#include <vector>
#include <chrono>
#include <unordered_map>
#include <fstream>

#include "boost/multi_array.hpp"
#include "boost/algorithm/string.hpp"
#include "batch_simulator.cpp"

int main(int argc, char** argv){
    if (argc != 2){
        cout << "Expected one argument corresponding to input tsv file." << endl;
        throw;
    }
    
    cout << "Opening file.." << endl;
    ifstream fin(argv[1]);

    string first_line;
    string line;
    vector<BatchSimulator> all_simulators;
    
    cout << "Checking first line for header.." << endl;

    // Check first line 
    //
    getline(fin, first_line);
    vector<string> check_first_line;
    split(check_first_line, first_line, boost::is_any_of("\t"));
    if (check_first_line[0] != "BatchID") {
        // Header not detected, so go back to start
        fin.seekg(0, ios::beg);
    }
    
    cout << "Parsing input file" << endl;

    while (getline(fin, line)) {
        // Split line into tab-separated parts
        BatchSimulator batch_simulator(line);
        batch_simulator.tsv_parser.Print();
        all_simulators.push_back(batch_simulator);
    }
    fin.close();
    
    string output_file_name = all_simulators[0].tsv_parser.output_file_name +".output.tsv";
    string summary_file_name = all_simulators[0].tsv_parser.output_file_name +".summary.tsv";

    vector<string> tsv_lines;
    vector<string> summary_lines;
    
    cout << "Simulating.." << endl;
    
    std::ofstream tsv_output;
    tsv_output.open(output_file_name);
    
    std::ofstream tsv_summary;
    tsv_summary.open(summary_file_name);

    tsv_output << TsvParser().OutputHeader() << endl;
    tsv_summary << TsvParser().SummaryHeader() << endl;

    for ( BatchSimulator simulator : all_simulators) {
        cout << "id: " << simulator.tsv_parser.data_id << endl;
        simulator.SimulateAll();
        tsv_lines.push_back(simulator.tsv_parser.StringifyOutput());
        tsv_output << simulator.tsv_parser.StringifyOutput() << endl;
        summary_lines.push_back(simulator.tsv_parser.StringifySummary());
        tsv_summary << simulator.tsv_parser.StringifySummary() << endl;
    }
    
    tsv_output.close();
    tsv_summary.close();

    return 0;
}

