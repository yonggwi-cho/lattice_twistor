#include <math>
#include "twistor.hpp"

#include "Eigen/Sparse"

Twistor::Twistor(int Ns):Ns(Ns),Nsp(2),Ndir(4){
  this->set_matrix();
  SparseMatrix<std::complex<double>> a;
}

void Twistor::set_matrix(void){
  int i,j;
  for(int alpha=0;alpha<this->Nsp;alpha++){
    for(int beta=0;beta<this->Nsp;beta++){
      for(int n1=0;n1<this->Ns;n1++){
	for(int n1p=0;n1p<this->Ns;n1p++){
	  for(int n2=0;n2<this->Ns;n2++){
	    for(int n2p=0;n2p<this->Ns;n2p++){
	      for(int n3=0;n3<this->Ns;n3++){
		for(int n3p=0;n3p<this->Ns;n3p++){
		  for(int n4=0;n4<this->Ns;n4++){
		    for(int n4p=0;n4p<this->Ns;n4p++){
		      i =  alpha + this->Nspn*n1 + this->Nspn*this->Ns*n2
			+ this->Nspn*pow(this->Ns*,2)*n3 + this->Nspn*pow(this->Ns,3)*n4;
		      j =  beta  + this->Nspn*n1p + this->Nspn*this->Ns*n2p
			+ this->Nspn*pow(this->Ns*,2)*n3p + this->Nspn*pow(this->Ns,3)*n4p;
		      d  = this->sgm(mu,alpha,beta);
		      d *= this->delivative(mu,n1,n1p,n2,n2p,n3,n3p,n4,n4p);
		      d *= this->jacobian(n1,n2,n3,n4);
		      this->Dirac(i,j) += d;
		    }
		  }
		}
	      }
	    }
	  }
	}   
      }
    }
  }
}

void Twistor::eigen_solver(){
  ComplexEigenSolver<MatrixXcf>solver(this->Dirac);
  cout<<solver.eigenvalues()<<endl;
}

complex<double> Twistor::distance(int n1, int n2, int n3, int n4){
  return double(pow(n1,2)+pow(n2,2)+pow(n3,2)+pow(n4,2)) ;
}

complex<double >Twistor::smg(int mu, int alpha, int beta){
  
}

