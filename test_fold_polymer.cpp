#include "fold_polymer_sim.cpp"
#include <eigen3/Eigen/Eigenvalues> // header file

using namespace std;

//* INT MAIN FOR DEBUGGING *//
void AllDebug() {
    int net_sz = 57;
    vector<int> initial_copy(net_sz,0);
        //{ 0, 0, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    vector<int> initial_template(net_sz,0);
        //{ 0, 0, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 2, 0, 0 };
    vector<int> initial_fold(net_sz,0);
        //{ 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
    for (int i = 3; i < net_sz-3; i++) {
        initial_template[i] = 1;
    }
    
    initial_copy[3] = 1;
    initial_fold[3] = 2;
    //for (int i = 3; i < net_sz/2; i++) {
    //    initial_copy[i] = 1;
    //}

    //for (int i = 3; i < net_sz/2; i++) {
    //   initial_fold[i] = 2;
    //}

    double G_pol_in = -0.7;
    double k_copy_in = 1;
    double k_fold_in = 50;
    double g_fold_in = 1;
    double g_u_in = 0;
    double g_h_in = 1;

    boostmat<double,2> energy_matrix_in;
    energy_matrix_in.resize(boost::extents[3][3]);

    energy_matrix_in[1][1] = 2;
    energy_matrix_in[2][1] = 0;
    energy_matrix_in[1][2] = 0;
    energy_matrix_in[2][2] = 2;

    int kinetic_flag_in = 0;
    
    
    double w = exp(g_h_in+g_fold_in);
    double v = exp(g_fold_in);
    double u = 1;

    Eigen::Matrix<double, 3, 3>  lifson_roig_transfer; // declare a real (double) 2x2 matrix
    lifson_roig_transfer << w, v, 0, 0, 0, u, v, v, u; // defined the matrix A

    Eigen::EigenSolver<Eigen::Matrix<double, 3,3> > eigen_lifson(lifson_roig_transfer); // the instance s(A) includes the eigensystem
            

    FoldDrivenSim fold_driven_sim(initial_copy, initial_template, initial_fold, G_pol_in - log(eigen_lifson.eigenvalues().real().maxCoeff()), k_copy_in,
        k_fold_in, g_fold_in, g_u_in, g_h_in, energy_matrix_in, kinetic_flag_in);

    random_device r{};
    unsigned int ran = r();
    for (int i = 0; i < 1000; i++) {
        cout << "Start " << i << endl;   
        fold_driven_sim.DebugSimulator(i);
        cout << "ACCURACY: " << fold_driven_sim.GetAccuracy() << endl;
    }
}

//* INT MAIN FOR DEBUGGING *//
void FoldOnlyDebug() {
    int net_sz = 26;
    vector<int> initial_copy(net_sz,0);
        //{ 0, 0, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    vector<int> initial_template(net_sz,0);
        //{ 0, 0, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 2, 0, 0 };
    vector<int> initial_fold(net_sz,0);
        //{ 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
    for (int i = 3; i < net_sz-3; i++) {
        initial_template[i] = 1;
    }

    for (int i = 3; i < net_sz/2; i++) {
        initial_copy[i] = 1;
    }

    for (int i = 3; i < net_sz/2; i++) {
        initial_fold[i] = 2;
    }

    double G_pol_in = 0.5;
    double k_copy_in = 0;
    double k_fold_in = 1;
    double g_fold_in = 0.8;
    double g_u_in = 0;
    double g_h_in = 0.1;

    boostmat<double,2> energy_matrix_in;
    energy_matrix_in.resize(boost::extents[3][3]);

    energy_matrix_in[1][1] = 1;
    energy_matrix_in[2][1] = 1;
    energy_matrix_in[1][2] = 1;
    energy_matrix_in[2][2] = 1;

    int kinetic_flag_in = 0;

    FoldDrivenSim fold_driven_sim(initial_copy, initial_template, initial_fold, G_pol_in, k_copy_in,
        k_fold_in, g_fold_in, g_u_in, g_h_in, energy_matrix_in, kinetic_flag_in);

    random_device r{};
    unsigned int ran = r();
    
    fold_driven_sim.DebugSimulator(1);
}

void FoldOnlyCoilFractions() {
    int net_sz = 2006;
    vector<int> initial_copy(net_sz,0);
        //{ 0, 0, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    vector<int> initial_template(net_sz,0);
        //{ 0, 0, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 2, 0, 0 };
    vector<int> initial_fold(net_sz,0);
        //{ 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
    for (int i = 3; i < net_sz-3; i++) {
        initial_template[i] = 1;
    }

    for (int i = 3; i < net_sz/2; i++) {
        initial_copy[i] = 1;
    }

    for (int i = 3; i < net_sz/2; i++) {
        initial_fold[i] = 1;
    }

    initial_fold[3] = 2;
    initial_fold[net_sz/2-1] = 2;

    double G_pol_in = 0.5;
    double k_copy_in = 0;
    double k_fold_in = 1;
    double g_fold_in = 0.1;
    double g_u_in = 1;
    double g_h_in = 0.8;

    boostmat<double,2> energy_matrix_in;
    energy_matrix_in.resize(boost::extents[3][3]);

    energy_matrix_in[1][1] = 1;
    energy_matrix_in[2][1] = 1;
    energy_matrix_in[1][2] = 1;
    energy_matrix_in[2][2] = 1;

    int kinetic_flag_in = 0;

    FoldDrivenSim fold_driven_sim(initial_copy, initial_template, initial_fold, G_pol_in, k_copy_in,
        k_fold_in, g_fold_in, g_u_in, g_h_in, energy_matrix_in, kinetic_flag_in);

    random_device r{};
    unsigned int ran = r();
    
    int J = 100;
    double avg_coil_num = 0;
    double avg_coil_tip = 0;
    for (int j = 0; j < J; j++) {
        cout << j << endl;
        vector<int> out_coil;
        vector<int> out_tip;
        fold_driven_sim.FoldCoilFraction( r(), 10, out_coil, out_tip);
        avg_coil_num = avg_coil_num + double(out_coil.back())/J;
        avg_coil_tip = avg_coil_tip + double(out_tip.back()-1)/J;
    }
    cout << (avg_coil_num-1)/(net_sz/2-4) << endl;
    cout << avg_coil_tip << endl;
}

void TestSpeed() {
    int net_sz = 507;
    vector<int> initial_copy(net_sz,0);
        //{ 0, 0, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    vector<int> initial_template(net_sz,0);
        //{ 0, 0, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 2, 0, 0 };
    vector<int> initial_fold(net_sz,0);
        //{ 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
    for (int i = 3; i < net_sz-3; i++) {
        initial_template[i] = 1;
    }
    
    initial_copy[3] = 1;
    initial_fold[3] = 2;
    //for (int i = 3; i < net_sz/2; i++) {
    //    initial_copy[i] = 1;
    //}

    //for (int i = 3; i < net_sz/2; i++) {
    //   initial_fold[i] = 2;
    //}

    double G_pol_in = -0.7;
    double k_copy_in = 1;
    double k_fold_in = 50;
    double g_fold_in = 1;
    double g_u_in = 0;
    double g_h_in = 1;

    boostmat<double,2> energy_matrix_in;
    energy_matrix_in.resize(boost::extents[3][3]);

    energy_matrix_in[1][1] = 2;
    energy_matrix_in[2][1] = 0;
    energy_matrix_in[1][2] = 0;
    energy_matrix_in[2][2] = 2;

    int kinetic_flag_in = 0;
    
    double w = exp(g_h_in+g_fold_in);
    double v = exp(g_fold_in);
    double u = 1;

    Eigen::Matrix<double, 3, 3>  lifson_roig_transfer; // declare a real (double) 2x2 matrix
    lifson_roig_transfer << w, v, 0, 0, 0, u, v, v, u; // defined the matrix A

    Eigen::EigenSolver<Eigen::Matrix<double, 3,3> > eigen_lifson(lifson_roig_transfer); // the instance s(A) includes the eigensystem
            
    cout << eigen_lifson.eigenvalues().real().maxCoeff() << endl;
    FoldDrivenSim fold_driven_sim(initial_copy, initial_template, initial_fold, G_pol_in - log(eigen_lifson.eigenvalues().real().maxCoeff()), k_copy_in,
        k_fold_in, g_fold_in, g_u_in, g_h_in, energy_matrix_in, kinetic_flag_in);
    
    random_device r{};
    unsigned int ran = r();
    for (int i = 0; i < 50; i++) {
        cout << i << endl;
        fold_driven_sim.Simulate(i);
        cout << fold_driven_sim.GetAccuracy() << endl;
    }
}

void TestShortAnalyticFolding() {
    int net_sz = 13;

    vector<int> initial_copy(net_sz,0);
        //{ 0, 0, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    vector<int> initial_template(net_sz,0);
        //{ 0, 0, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 2, 0, 0 };
    vector<int> initial_fold(net_sz,0);
        //{ 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
    for (int i = 3; i < net_sz-3; i++) {
        initial_template[i] = 1;
    }

    initial_fold[3] = 2;
    initial_copy[3] = 1;

    double G_pol_in = 0;
    double k_copy_in = 1;
    double k_fold_in = 0;
    double g_fold_in = 0;
    double g_u_in = 0;
    double g_h_in = 0;

    boostmat<double,2> energy_matrix_in;
    energy_matrix_in.resize(boost::extents[3][3]);

    energy_matrix_in[1][1] = 2;
    energy_matrix_in[2][1] = 0;
    energy_matrix_in[1][2] = 0;
    energy_matrix_in[2][2] = 2;

    int kinetic_flag_in = 0;

    FoldDrivenSim fold_driven_sim(initial_copy, initial_template, initial_fold, G_pol_in, k_copy_in,
        k_fold_in, g_fold_in, g_u_in, g_h_in, energy_matrix_in, kinetic_flag_in);

    random_device r{};
    unsigned int ran = r();
    
    int J = 100000;
    int count = 0;
    double avg_coil_num = 0;
    double avg_coil_tip = 0;
    for (int j = 0; j < J; j++) {
        fold_driven_sim.Simulate( r());
        if ( (fold_driven_sim.complete_copy[4] == 1) + (fold_driven_sim.complete_copy[5] == 1) + (fold_driven_sim.complete_copy[6] == 1) 
               + (fold_driven_sim.complete_copy[7] == 1) + (fold_driven_sim.complete_copy[8] == 1) + (fold_driven_sim.complete_copy[9] == 1) == 6) {
            count = count + 1;
        }
    }
    cout << count << endl;
}

int main() {
    TestShortAnalyticFolding();
}