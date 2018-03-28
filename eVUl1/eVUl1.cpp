// eVUl1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
using namespace Eigen;
double qfun(const Ref<const VectorXd>& x, const Ref<const VectorXd>& b, const Ref<const MatrixXd> &A);
void qGradient(const Ref<const VectorXd>& x, const Ref<const VectorXd>& b, const Ref<const MatrixXd> &A,
	Ref<VectorXd> s);
struct Param
{
	double eMax;
	double sigma;
	double epsilon;
};

int main()
{
	std::string line;
	double d1, mu, tau, alpha, eps, t0, t3, t4, t5;
	double t1 = 0;
	double t2 = 0;
	mu = 4;
	int n1;
	char* end;
	char* e1;
	Param pr = { 1,0.01,1e-6 };
	std::size_t i;
	int j = -1;
	std::size_t rows = 862;
	std::size_t cols = 2;
	MatrixXd A(rows, cols);
	SparseMatrix<double> U(cols, cols);//U is an n-by-n matrix for storage of the basis matrix U_k with dimension n_k
	U.reserve(VectorXi::Constant(cols, 1));
	VectorXd b(rows);
	VectorXd x(cols);
	VectorXd xplus(cols);
	VectorXd w(cols);//used in the shrink operator
	VectorXd dqx(cols);//gradient of q(x) at x
	VectorXd dqp(cols);//gradient of q(x) at p
	VectorXd g(cols);//g^k in the algorithm
	VectorXd ud(cols);// u^*_k in the algorithm
	MatrixXd Q;//Q_k inverse in the algorithm
	VectorXd p(cols);//p^k in the algorithm
	p.setZero();
	x.setOnes();//initial value x0
	ud.setZero();
	//Read input data with SVMlight format
	std::fstream input("../Datasets/fourclass.txt", std::ios::in);
	for (i = 0; i < rows; i++)
	{
		std::getline(input, line);
		b(i) = strtod(line.c_str(), &end);
		while (1)
		{
			n1 = strtoul(end, &e1, 10);
			d1 = strtod(e1 + 1, &end);
			if (d1 == 0)
				break;
			A(i, n1 - 1) = d1;                    // alternative for sparse A: mat.coeffRef(i,j) += v_ij;
		}
	}
	qGradient(x, b, A, dqx);
	VectorXd vft = A.transpose()*b;
	tau = 0.1* vft.lpNorm<Infinity>();
	//====Shrink operator plus calculating epsilon beginning======
	do {
		alpha = tau / mu;
		w = x - 1 / mu * dqx;
		eps = tau * x.lpNorm<1>();// initialize eps
		for (i = 0; i < cols; i++)
		{
			if (w(i) > alpha)
			{
				t1 -= x(i);
				p(i) = w(i) - alpha;
			}
			else
				if (w(i) < -alpha)
				{
					t1 += x(i);
					p(i) = w(i) + alpha;
				}
				else
					t2 += x(i)*(dqx(i) - mu * x(i));
		}
		t1 *= tau;
		eps = eps + t1 + t2;//the value of eps for usage
		if (eps < 0)
			std::cerr << "=================eps=%d <0! Wrong!=================" << eps << std::endl;
		//====Shrink operator  end======
		qGradient(p, b, A, dqp);
		g = mu * (x - p) + dqp - dqx;

		// Calculating U_k
		while (j < 0)
		{
			for (i = 0; i < cols; i++)
				if (abs(x(i)) > eps / 2)
				{
					j++;
					U.insert(i, j) = 1;
				}
			eps /= 3;// Consider moving this line 
		}
		// U-step
		Q = U.leftCols(j + 1).transpose() *A.transpose()*A * U.leftCols(j + 1);
		ud.segment(0, j + 1) = Q.llt().solve(-U.leftCols(j + 1).transpose() * g);
		xplus = p + U.leftCols(j + 1)*ud.segment(0, j + 1);
		// Serious step test
		t0 = tau * (x.lpNorm<1>() - xplus.lpNorm<1>());
		t3 = qfun(x, b, A) - qfun(xplus, b, A) + t0;
		t4 = dqx.transpose()*(x - p) + t0;
		if (t3 >= pr.sigma*t4)
		{
			x = xplus;
			qGradient(x, b, A, dqx);
			std::cout << "serious step!\n";
			// Later, we need to update mu
			// mu= 
		}
		else
			mu *= 2;
		//Update These values are going to be reevaluated
		j = -1;
		t1 = 0;
		t2 = 0;
		U.setZero();// Later maybe think about not recalculating U_k from the beginning
		ud.setZero();
		p.setZero();
	} while (g.norm() > pr.epsilon);// Later: this stopping criterion is not good
	system("pause");// for debugging
	return 0;
}

void qGradient(const Ref<const VectorXd>& x, const Ref<const VectorXd>& b, const Ref<const MatrixXd> &A,
	Ref<VectorXd> s)
{
	s = A.transpose() *(A*x - b);

}
double qfun(const Ref<const VectorXd>& x, const Ref<const VectorXd>& b, const Ref<const MatrixXd> &A)
{
	double qf = (x.transpose()*A.transpose() - b.transpose())*(A*x - b);
	return qf / 2;

}
