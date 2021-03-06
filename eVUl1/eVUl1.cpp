#include "stdafx.h"
using namespace Eigen;
using namespace std;
namespace fs = std::experimental::filesystem::v1;
double qfun(const Ref<const VectorXd>& x
	//, const Ref<const VectorXd>& b, const Ref<const MatrixXd> &A
);
void qGradient(const Ref<const VectorXd>& x
	//, const Ref<const VectorXd>& b, const Ref<const MatrixXd> &A
	, Ref<VectorXd> s);
void Dh(const Ref<const VectorXd>& x, Ref<VectorXd> s);
struct Param
{
	double tau;
	double sigma;
	double epsilon;
	double mu;
};
double eVUl1(Ref<VectorXd>& x, const Ref<const MatrixXd> &M,
	double(*qfun)(const Ref<const VectorXd>& x,
		//const Ref<const VectorXd>& b, const Ref<const MatrixXd> &A),
		void(*qGradient)(const Ref<const VectorXd>& x,
			//const Ref<const VectorXd>& b, const Ref<const MatrixXd> &A,
			Ref<VectorXd> s), Param pr));//later consider using sparse matrix M
VectorXd x;
MatrixXd A;
VectorXd b;
VectorXd w;//used in the shrink operator
VectorXd dqx;//gradient of q(x) at x
VectorXd dqp;//gradient of q(x) at p
VectorXd g;//g^k in the algorithm
VectorXd ud;// u^*_k in the algorithm
VectorXd p;//p^k in the algorithm
VectorXd xplus;
MatrixXd Q;//Q_k inverse in the algorithm
VectorXd xkm = x;//used in the update of mu 
VectorXd dhxkm;
VectorXd dqxkm;//used in the update of mu 
VectorXd gxkm;//used in the update of mu 
VectorXd t1;//used in the update of mu 
VectorXd dhx;

int main()
{	//string line;
	//double d1, mu, tau, alpha, eps, t0, t3, t4, t5, mx_mp;//mx_mp is m(x)-m(p) in the algorithm
	//double t1 = 0;
	//double t2 = 0;
	//mu = 4;
	//int n1, k;
	//char* end;
	//char* e1;
	//k = 0;//iteration index
	//Param pr = { 1,0.01,1e-6 };
	//size_t i;
	//int j = -1;
	//size_t rows = 3;
	//size_t n = 3;
	//MatrixXd A(rows, n);
	//SparseMatrix<double> U(n, n);//U is an n-by-n matrix for storage of the basis matrix U_k with dimension n_k
	//U.reserve(VectorXi::Constant(n, 1));
	//VectorXd b(rows);
	//A << 2, 0, 0,
	//	0, 1, 0,
	//	0, 0, 0.01;
	//b << 1, 0.5, 0.1;
	//VectorXd x(n);
	//VectorXd xplus(n);
	//VectorXd w(n);//used in the shrink operator
	//VectorXd dqx(n);//gradient of q(x) at x
	//VectorXd dqp(n);//gradient of q(x) at p
	//VectorXd g(n);//g^k in the algorithm
	//VectorXd ud(n);// u^*_k in the algorithm
	//MatrixXd Q;//Q_k inverse in the algorithm
	//VectorXd p(n);//p^k in the algorithm
	//p.setZero();
	//x.setOnes();//initial value x0
	//x *= 0.1;
	//ud.setZero();
	////Read input data with SVMlight format
	//// fstream input("../Datasets/bodyfat",  ios::in);
	////for (i = 0; i < rows; i++)
	////{
	////	 getline(input, line);
	////	b(i) = strtod(line.c_str(), &end);
	////	// cout << b(i);
	////	while (1)
	////	{
	////		n1 = strtoul(end, &e1, 10);
	////		// cout << n1<<endl;
	////		if (n1 == 0)
	////			break;
	////		d1 = strtod(e1 + 1, &end);
	////		// cout << d1<<endl;			
	////		A(i, n1 - 1) = d1;                    // alternative for sparse A: mat.coeffRef(i,j) += v_ij;
	////	}
	////}
	//// cout << b<<endl;
	//MatrixXd M = A.transpose()*A;
	//qGradient(x, b, A, dqx);
	//VectorXd vft = A.transpose()*b;
	//tau = 0.1* vft.lpNorm<Infinity>();
	////====Shrink operator plus calculating epsilon beginning======
	//do {
	//	alpha = tau / mu;
	//	w = x - 1 / mu * dqx;
	//	//eps = tau * x.lpNorm<1>();// initialize eps
	//	for (i = 0; i < n; i++)
	//	{
	//		if (w(i) > alpha)
	//		{
	//			//t1 -= x(i);
	//			p(i) = w(i) - alpha;
	//		}
	//		else
	//			if (w(i) < -alpha)
	//			{
	//				//t1 += x(i);
	//				p(i) = w(i) + alpha;
	//			}
	//		//else
	//			//t2 += x(i)*(dqx(i) - mu * x(i));
	//	}
	//	// cout << "p="<<p << endl;
	//	mx_mp = dqx.transpose()*(x - p) + tau * (x.lpNorm<1>() - p.lpNorm<1>());
	//	eps = mx_mp - mu * (p - x).transpose()*(p - x);
	//	//double re = 0.4 + 0.2*0.2 - (0.01 + pow(0.35,2 ))/2-0.2*0.6-4*(pow(0.35, 2) +pow(0.05,2));
	//	//t1 *= tau;
	//	// cout << "re=" << re << endl;
	//	//eps = eps + t1 + t2;//the value of eps for usage
	//	cout << "eps original=" << eps << endl;
	//	if (eps < 0)
	//		cerr << "=================eps=%d <0! Wrong!=================" << eps << endl;
	//	//====Shrink operator  end======
	//	qGradient(p, b, A, dqp);
	//	g = mu * (x - p) + dqp - dqx;
	//	// Calculating U_k
	//	for (i = 0; i < n; i++)
	//		if (abs(x(i)) > eps / 2)
	//		{
	//			j++;
	//			U.insert(i, j) = 1;
	//		}
	//	if (j < 0)
	//	{
	//		xplus = p;
	//	}
	//	else
	//	{
	//		//eps /= 3;// Consider moving this line 
	////}
	////cout <<"eps="<< eps<<"and j ="<<j << endl;
	////cout << U.leftCols(j + 1);
	//// U-step
	//		Q = U.leftCols(j + 1).transpose() *M* U.leftCols(j + 1);
	//		ud.segment(0, j + 1) = Q.llt().solve(-U.leftCols(j + 1).transpose() * g);
	//		xplus = p + U.leftCols(j + 1)*ud.segment(0, j + 1);
	//	}
	//	// Serious step test
	//	t0 = tau * (x.lpNorm<1>() - xplus.lpNorm<1>());
	//	t3 = qfun(x, b, A) - qfun(xplus, b, A) + t0;
	//	t5 = t3 / mx_mp;
	//	if (isnan(t5))
	//	{
	//		cout << "===============t5 is not a number!=================\n";
	//		VectorXd v4 = x - p;
	//		cout << "x=" << x << '\n';
	//		cout << "p=" << p << '\n';
	//		cout << v4 << '\n';
	//	}
	//	if (t3 >= pr.sigma*mx_mp)
	//	{
	//		x = xplus;
	//		qGradient(x, b, A, dqx);
	//		cout << "serious step!\n";
	//		// Later, we need to update mu
	//		// mu= 
	//	}
	//	else
	//		mu *= 2;
	//	cout << "Iteration " << k << ": dim Uk -1 = " << j << endl;
	//	//Update These values are going to be reevaluated
	//	j = -1;
	//	t1 = 0;
	//	t2 = 0;
	//	U.setZero();// Later maybe think about not recalculating U_k from the beginning
	//	ud.setZero();
	//	p.setZero();
	//	k++;
	//} while (g.norm() > pr.epsilon);// Later: this stopping criterion is not good
	//cout << "returned value is f(x) or f(x_+), and f(x)=" << qfun(x, b, A) + tau * x.lpNorm<1>() << '\n';
	//cout << "f(x_+)=" << qfun(xplus, b, A) + tau * xplus.lpNorm<1>() << '\n';
	//// cout << "optimal value is norm(b)^2 /2=" << b.transpose()*b / 2 << '\n';

//fs::path dataPath = "C:/Users/Leo/Source/eVUl1/Datasets";
//std::string s;
//std::cout << "Maximum size of a string is " << s.max_size() << "\n";
	fs::path dataPath = "../Datasets/";
	//cout << p2;
	unsigned long i, n1;
	double d1;
	string line;
	char* end;
	char* e1;
	for (auto& p : fs::directory_iterator(dataPath))
	{
		std::cout << p << '\n';
		////Read input data with SVMlight format
		fstream input(p, ios::in);
		i = 0;
		for (line; getline(input, line); )
		{
			b(i) = stod(line);
			// cout << b(i);
			while (1)
			{
				n1 = stoul(line);
				// cout << n1<<endl;
				if (n1 == 0)
					break;
				d1 = stod(line);
				// cout << d1<<endl;			
				A(i, n1 - 1) = d1;                    // alternative for sparse A: mat.coeffRef(i,j) += v_ij;
			}
			i++;
		}
		cout << b;

	}
	system("pause");// for debugging
	return 0;
}

void qGradient(const Ref<const VectorXd>& x,
	//const Ref<const VectorXd>& b, const Ref<const MatrixXd> &A,
	Ref<VectorXd> s)
{
	s = A.transpose() *(A*x - b);

}
double qfun(const Ref<const VectorXd>& x
	//, const Ref<const VectorXd>& b, const Ref<const MatrixXd> &A
)
{
	double qf = (x.transpose()*A.transpose() - b.transpose())*(A*x - b);
	return qf / 2;

}
void Dh(const Ref<const VectorXd>& x, Ref<VectorXd> s) {
	unsigned long i;
	size_t n = x.size();
	s.setZero();
	for (i = 0; i < n; i++)
		if (x(i) > 0)
			s(i) = 1;
		else
			if (x(i) < 0)
				s(i) = -1;
}
double eVUl1(Ref<VectorXd>& x, const Ref<const MatrixXd> &M,
	double(*qfun)(const Ref<const VectorXd>& x),
	void(*qGradient)(const Ref<const VectorXd>& x, Ref<VectorXd> s), Param pr)
{
	unsigned long i, k = 0, j = 0;
	size_t n = x.size();
	double alpha, mu = pr.mu, mx_mp, eps, x1norm, fx, fxplus, fv, t, mun;
	//VectorXd w(n);//used in the shrink operator
	//VectorXd dqx(n);//gradient of q(x) at x
	//VectorXd dqp(n);//gradient of q(x) at p
	//VectorXd g(n);//g^k in the algorithm
	//VectorXd ud(n);// u^*_k in the algorithm
	//VectorXd p(n);//p^k in the algorithm
	//VectorXd xplus(n);
	//MatrixXd Q;//Q_k inverse in the algorithm
	//VectorXd dqxkm;//used in the update of mu 
	//VectorXd gxkm;//used in the update of mu 
	//VectorXd t1;//used in the update of mu 
	//U is an n-by-n matrix for storage of the basis matrix U_k with dimension n_k 
	SparseMatrix<double> U(n, n);
	U.reserve(VectorXi::Constant(n, 1));
	p.setZero();
	x1norm = x.lpNorm<1>();
	fx = qfun(x) + pr.tau*x1norm;
	qGradient(x, dqx);
	xkm = x;//used in the update of mu 
	Dh(x, dhx);
	dqxkm = dqx;
	dhxkm = dhx;

	do {
		alpha = pr.tau / mu;
		w = x - 1 / mu * dqx;
		p = (w.array() > alpha).select(w.array() - alpha, p);
		p = (w.array() < -alpha).select(w.array() + alpha, p);
		mx_mp = dqx.transpose()*(x - p) + pr.tau*(x1norm - p.lpNorm<1>());
		eps = mx_mp - mu * (p - x).transpose()*(p - x);
		cout << endl << "Iter. " << k << ": epsilon = " << eps;
		//fprintf('\nIter. %d: epsilon =%e; ', k, eps);
		if (eps < 0)
			cout << endl << "Warning: epsilon <0, something is wrong!" << endl;
		//fprintf('\nWarning: \n');
		qGradient(p, dqx);
		g = mu * (x - p) + dqp - dqx;
		// Calculating U_k
		for (i = 0; i < n; i++)
			if (abs(x(i)) > eps / 2)
			{
				U.insert(i, j) = 1;
				j++;
			}
		if (j < 1)
		{
			xplus = p;
		}
		else
		{
			Q = U.leftCols(j).transpose() *M* U.leftCols(j);
			ud.segment(0, j) = Q.llt().solve(-U.leftCols(j).transpose() * g);
			xplus = p + U.leftCols(j)*ud.segment(0, j);

		}
		//('%d \n', j);

		fxplus = qfun(xplus) + pr.tau*xplus.lpNorm<1>();
		if (fx - fxplus >= pr.sigma*mx_mp)
		{// Update mu
			t = max(0.04, 0.1*mu);
			gxkm = dqxkm + pr.tau * dhxkm;
			t1 = dqx + pr.tau * dhx - gxkm;
			mun = t1.transpose()*t1;
			mun = mun / (t1.transpose()*(x - xkm));
			t = max(mun, t);
			mu = min(10 * mu, t);
			xkm = x;
			dqxkm = dqx;
			dhxkm = dhx;
			x = xplus;
			x1norm = x.lpNorm<1>();
			fx = qfun(x) + pr.tau*x1norm;
			qGradient(x, dqx);
			k++;
		}
		else
			mu *= 2;
		cout << " dim Uk=" << j << endl;
		j = 0;
		U.setZero();
		p.setZero();
		ud.setZero();
	} while (g.norm() > pr.epsilon);
	fv = min(fx, fxplus);
	cout << "eVUl1 is terminated successfully with final value of the objective function: " << fv << endl;
	return fv;
}