
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

const int dim=2;
void step(double *p, double *q, double dt, const int dim, double & H);
void init (double *p, double *q, double & H, const double e);
int main()
{
	const double PI=3.1415926536;
	double dt=0.05, t=0., t_end=20*PI;
	//const int dim=2;
	double p[dim], q[dim];
	const double e=0.6;
	double H;
	const int N=t_end/dt; //const int N=125000;
	
	//init conditions
	ofstream out("symp_euler1.txt");
	init (p, q,  H,  e);
	out<<t<<"\t "<< q[0]<<"\t "<< q[1]<<"\t "<< p[0]<<"\t "<< p[1]<<"\t "<< H<<endl;
	//out<<q[0]<<"\t "<<q[1]<<endl;

	//next step

	for(int i=0; i<N; i++)
	{
		
		step(p, q, dt,  dim,   H);
		t+=dt;
		out<<t<<"\t"<<q[0]<<"\t"<<q[1]<<"\t"<<p[0]<<"\t"<<p[1]<<"\t"<<H<<endl;
		//out<<q[0]<<"\t   "<<   q[1]<<endl;
	}
	out.close();
  
    return 0;
}

void step(double *p, double *q, double dt, const int dim, double &  H)
{
	double q2=q[0]*q[0]+q[1]*q[1];
	for(int i=0; i<dim; i++)
	{
		p[i]=p[i]-dt*q[i]/pow(q2,3.0/2.0);
		q[i]=q[i]+dt*p[i];
	}
	H=0.5*(p[0]*p[0]+p[1]*p[1])-1.0/sqrt(q[0]*q[0]+q[1]*q[1]);
}

void init (double *p, double *q, double & H, const double e)
{
	p[0]=0;	p[1]=sqrt((1+e)/(1-e));
	q[0]=1-e;	q[1]=0;
	
	H=0.5*(p[0]*p[0]+p[1]*p[1])-1.0/sqrt(q[0]*q[0]+q[1]*q[1]);
}
