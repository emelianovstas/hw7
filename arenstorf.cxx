#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

void control(double* y, double* rk4, double* rk5, const double tol, double& t);
void func(double* f, const double* const y, const double _mu)
{
double r=sqrt((y[0]+_mu)*(y[0]+_mu)+y[2]*y[2]);
double s=sqrt((y[0]-1+_mu)*(y[0]-1+_mu)+y[2]*y[2]);

	f[0]=y[1];
	f[2]=y[3];
	f[1]=y[0]+2*y[3]-((1-_mu)*(y[3]+_mu))/(pow(r,3))-(_mu*(y[0]-1+_mu))/(pow(s,3));
	f[3]=y[2]-2*y[1]-((1-_mu)*y[2])/(pow(r,3))-(_mu*y[2])/(pow(s,3));
}
int main()
{
	const int dim=4;
	double dt=1e-3, dtn, t=0;
	const double T=17.056216560157;
	const double _mu=0.012277571;
	const double tol=1e-3;
	double rkstep[dim];
	double rkmax;

	double y[dim];
	y[0]=0.994;
	y[1]=0;
	y[2]=0;
	y[3]=-2.0015810637908;
	
	double rk4[dim];
	double rk5[dim];
	double k1[dim], k2[dim], k3[dim], k4[dim], k5[dim], k6[dim], k7[dim];
	double ytemp[dim];
// ytemp1[dim],ytemp2[dim],ytemp3[dim],ytemp4[dim],ytemp5[dim],ytemp6[dim];
int i=0;
ofstream out("out.txt");
out<<t<<"\t"<<dt<<"\t"<<y[0]<<"\t"<<y[1]<<endl;
while(T>t)
{
i++;
//cout<<t<<"\t"<<dt<<"\t"<<y[0]<<"\t"<<y[1]<<endl;

t+=dt;
//cout<<t<<"\t"<<dt<<"\t"<<y[0]<<"\t"<<y[1]<<endl;

func(k1,y,_mu);
	ytemp[0] = y[0] + dt/5.*k1[0];
	ytemp[1] = y[1] + dt/5.*k1[1];
	ytemp[2] = y[2] + dt/5.*k1[2];
	ytemp[3] = y[3] + dt/5.*k1[3];
func(k2,ytemp,_mu);
	ytemp[0] = y[0] + dt*(3./40.*k1[0] + 9./40.*k2[0]);
	ytemp[1] = y[1] + dt*(3./40.*k1[1] + 9./40.*k2[1]);
	ytemp[2] = y[2] + dt*(3./40.*k1[2] + 9./40.*k2[2]);
	ytemp[3] = y[3] + dt*(3./40.*k1[3] + 9./40.*k2[3]);
func(k3,ytemp,_mu);
	ytemp[0] = y[0] + dt * (44./45.*k1[0] - 56./15.*k2[0] + 32./9.*k3[0]);
	ytemp[1] = y[1] + dt * (44./45.*k1[1] - 56./15.*k2[1] + 32./9.*k3[1]);
	ytemp[2] = y[2] + dt * (44./45.*k1[2] - 56./15.*k2[2] + 32./9.*k3[2]);
	ytemp[3] = y[3] + dt * (44./45.*k1[3] - 56./15.*k2[3] + 32./9.*k3[3]);
func(k4,ytemp,_mu);
	ytemp[0] = y[0] + dt * (19372./6561.*k1[0] - 25360./2187.*k2[0] + 64448./6561.*k3[0] - 212./729.*k4[0]);
	ytemp[1] = y[1] + dt * (19372./6561.*k1[1] - 25360./2187.*k2[1] + 64448./6561.*k3[1] - 212./729.*k4[1]);
	ytemp[2] = y[2] + dt * (19372./6561.*k1[2] - 25360./2187.*k2[2] + 64448./6561.*k3[2] - 212./729.*k4[2]);
	ytemp[3] = y[3] + dt * (19372./6561.*k1[3] - 25360./2187.*k2[3] + 64448./6561.*k3[3] - 212./729.*k4[3]);
func(k5,ytemp,_mu);
	ytemp[0] = y[0] + dt * (9017./3168.*k1[0] - 355./33.*k2[0] + 46732./5247.*k3[0] + 49./176.*k4[0] - 5103./18656.*k5[0]);
	ytemp[1] = y[1] + dt * (9017./3168.*k1[1] - 355./33.*k2[1] + 46732./5247.*k3[1] + 49./176.*k4[1] - 5103./18656.*k5[1]);
	ytemp[2] = y[2] + dt * (9017./3168.*k1[2] - 355./33.*k2[2] + 46732./5247.*k3[2] + 49./176.*k4[2] - 5103./18656.*k5[2]);
	ytemp[3] = y[3] + dt * (9017./3168.*k1[3] - 355./33.*k2[3] + 46732./5247.*k3[3] + 49./176.*k4[3] - 5103./18656.*k5[3]);
func(k6,ytemp,_mu);
	ytemp[0] = y[0] + dt * (35./384.*k1[0] + 500./113.*k3[0] + 125./192.*k4[0] - 2187./6784.*k5[0] + 11./84.*k6[0]);
	ytemp[1] = y[1] + dt * (35./384.*k1[1] + 500./113.*k3[1] + 125./192.*k4[1] - 2187./6784.*k5[1] + 11./84.*k6[1]);
	ytemp[2] = y[2] + dt * (35./384.*k1[2] + 500./113.*k3[2] + 125./192.*k4[2] - 2187./6784.*k5[2] + 11./84.*k6[2]);
	ytemp[3] = y[3] + dt * (35./384.*k1[3] + 500./113.*k3[3] + 125./192.*k4[3] - 2187./6784.*k5[3] + 11./84.*k6[3]);
func(k7,ytemp,_mu);
//cout<<t<<"\t"<<dt<<"\t"<<y[0]<<"\t"<<y[1]<<endl;

for (int a=0; a<dim; a++)

	rk5[a] = y[a] + dt * (35./384.*k1[a] + 500./1113.*k3[a] + 125./192.*k4[a] - 2187./6784.*k5[a] + 11./48.*k6[a]);

for (int a=0;a<dim;a++)

	rk4[a] = y[a] + dt * (5179./57600.*k1[a] + 7571./16695.*k3[a] + 393./640.*k4[a] - 92097./339200.*k5[a] + 187./2100.*k6[a] + 1./40.*k7[a]); 

//control(y,rk5,tol,dt);
//t+=dt;
//cout<<t<<"\t"<<dt<<"\t"<<y[0]<<"\t"<<y[1]<<endl;

for (int a=0;a<dim;a++)
{
	rkstep[a]=abs(rk4[a]-rk5[a]);
}
	rkmax=rkstep[0];
for (int a=0;a<dim;a++)
{
		if (rkstep[a]>rkmax)
			rkmax=rkstep[a];
			else continue;
}

	dtn=dt*pow((tol/rkmax),1./5.)*0.8;	
	dt=dtn;

for(int a=0;a<dim;a++)

	y[a]=rk4[a];


out<<t<<"\t"<<dt<<"\t"<<y[0]<<"\t"<<y[1]<<"\t"<<y[2]<<"\t"<<y[3]<<endl;

}
out.close();
return 0;
}
 
