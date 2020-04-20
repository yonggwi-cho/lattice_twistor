#include<stdio>
#include<iostream>

#include "params.hpp"
#include "Eigen/Sparse"

using namespace std;
using namespace Eigen;

class Twistor {
public:
  Twistor(int Ns);
  void set_matrix();
  void eigen_solver(Matrix a);
private:
  const int Nsp,Ndir,Ns,Nvol,Nint;
  complex<double> *eigen_value;
  complex<double> *eigen_state;
  Matrix<complex<double>,this->Ns,this->Ns> Dirac;
  comple<double> Jacobian();
  comple<double> delivative();
  comple<double> distance();
}
