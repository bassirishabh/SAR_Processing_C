#include <iostream>
#include "Eigen/Dense"
#include <vector>
#include <math.h> 

using namespace std;

int main()
{
   // Eigen::MatrixXf A = Eigen::MatrixXf::Random(3, 2);
   // std::cout << "A: " << A << std::endl;

   // Eigen::JacobiSVD<Eigen::MatrixXf> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
   // std::cout << "U: " << svd.matrixU() << std::endl;;
   // std::cout << "V: " << svd.matrixV() << std::endl;;
   // std::cout << "Sigma: " << svd.singularValues() << std::endl;;
   int M=5;
   int N=10;
   int K= 5;

   Eigen::MatrixXd matched_filter;   // yet to take conjugate and complex numbers.
   matched_filter.setRandom(420,1);
   std::cout << "Here is the matrix m:\n" << matched_filter << std::endl;

   Eigen::MatrixXd time_delay_arr_P;
   time_delay_arr_P.setRandom(N,M);

   Eigen::MatrixXd slow_time_idx_arr_P;
   slow_time_idx_arr_P.setRandom(M,N);

   Eigen::MatrixXd range_cell_arr_P;
   range_cell_arr_P.setRandom(N,M);

   Eigen::MatrixXd Pnu;
   Pnu.setRandom(N*K,1);

   Eigen::MatrixXd nu;
   nu.setRandom(M,1);

   Eigen::MatrixXd ph_cor_coeffc;
   ph_cor_coeffc.setRandom(N,1);

   Eigen::MatrixXd result_matrix;
   result_matrix.setRandom(420,1);

   // Eigen::MatrixXd idx_values::Zero(420,N);
   // // idx_values.setRandom(420,N);
   // cout<< idx_values;

   Eigen::ArrayXXf idx_values = Eigen::ArrayXXf::Zero(420, N);

//    std::complex<double> c;
//    Eigen::MatrixXd R1; R1.setRandom(2,2);
//    Eigen::MatrixXcd C1 = c*R1; // multiply complex*real
//    Eigen::MatrixXcd C2 = c*C1; // complex scalar times complex matrix
//    C1(0,0) = c; // assign complex value.

   for(int m=0;m<M;m++){
         // std::complex<double> ph_cor_coeffc = exp (c * 2 * 3 * 1 * m);
         // cout<<"k"<<time_delay_arr_P(Eigen::all,m);
         for(int r=0;r<N;r++){
            // cout<<Pnu(idx_values(Eigen::all,r));
            Pnu(idx_values(Eigen::all,r)) += ph_cor_coeffc(r) * matched_filter * nu(m);
            
            // result_matrix+=ph_cor_coeffc(r) * matched_filter * nu(m);
            // for(int l=0;l<420;l++){
            //     int s=idx_values(l,r);
            //     Pnu(s) += ph_cor_coeffc(r) * matched_filter(l) * nu(m);
            // }
            
            // cout<<"Hello";
         }
   }
}
