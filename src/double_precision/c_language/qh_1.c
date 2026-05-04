#include <stdio.h>
#include <math.h>

double Qh01(double f(double),double h){
  double x, Ex, Epx, fx, S, SS;
  int n;
  x=0;
  Ex=exp(-exp(x*(1+exp(-x))));
  Epx=Ex*exp(x*(1+exp(-x)))*(1+(1-x)*exp(-x));
  S=h*f(Ex)*Epx;
  n=-1; x=-h;
  Ex=exp(-exp(x*(1+exp(-x))));
  Epx=Ex*exp(x*(1+exp(-x)))*(1+(1-x)*exp(-x));
  fx=f(Ex)*Epx;
  SS=S+h*fx;
  while((SS!=S)||(n>-10)){
    n-=1; x-=h;
    S=SS;
    Ex=exp(-exp(x*(1+exp(-x))));
    Epx=Ex*exp(x*(1+exp(-x)))*(1+(1-x)*exp(-x));
    fx=f(Ex)*Epx;
    SS +=h*fx;
  }
  n=1; x=h; 
  Ex=exp(-exp(x*(1+exp(-x))));
  Epx=Ex*exp(x*(1+exp(-x)))*(1+(1-x)*exp(-x));
  fx=f(Ex)*Epx;
  SS=S+h*fx;
  while((SS!=S)||(n<10)){
    n++; x+=h;
    S=SS;
    Ex=exp(-exp(x*(1+exp(-x))));
    Epx=Ex*exp(x*(1+exp(-x)))*(1+(1-x)*exp(-x));
    fx=f(Ex)*Epx;
    SS +=h*fx;
  }
  return S;
}

double Qhinfinf(double f(double),double h){
  double g(double x){
    double x1mx=x*(1-x);
    double x1mx2=x1mx*x1mx;
    return f((x-0.5)/x1mx)*(0.5-x1mx)/x1mx2;
  }
  return Qh01(g,h);
}

double Qhainf(double f(double), double a, double h){
  double g(double x){
    return f(a-1+1/x)/x/x;
  }
  return Qh01(g,h);
}

double Qhinfb(double f(double), double b, double h){
  double g(double x){
    return f(b+1-1/x)/x/x;
  }
  return Qh01(g,h);
}

double f1(double x){
  return 1;
}

double f2(double x){
  return x;
}

double f3(double x){
  return pow(x,99);
}

double f4(double x){
  return exp(x);
}

double f5(double x){
  return 1/(1+x*x) ;
}

double f6(double x){
  return exp(-M_PI*x*x);
}

int main(void){
  printf("1, f1, %.15e\n",Qh01(f1,pow(2,-6))-1);
  printf("2, f2, %.15e\n",Qh01(f2,pow(2,-6))*2-1);
  printf("3, f3, %.15e\n",Qh01(f3,pow(2,-6))*100-1);
  printf("4, f4, %.15e\n",Qh01(f4,pow(2,-6))/(exp(1)-1)-1);
  printf("5, f5, %.15e\n",Qhinfinf(f5,pow(2,-6))/M_PI-1);
  printf("6, f6, %.15e\n",Qhinfinf(f6,pow(2,-6))-1);
  printf("7, f5, %.15e\n",Qhainf(f5,0,pow(2,-6))/M_PI*2-1);
  printf("8, f5, %.15e\n",Qhainf(f5,1,pow(2,-6))/M_PI*4-1);
  printf("9, f6, %.15e\n",Qhinfb(f6,0,pow(2,-6))-0.5);
  printf("10, f6, %.15e\n",Qhinfb(f6,1,pow(2,-6))/((1+erf(sqrt(M_PI)))/2)-1);
  return 0;
}

/* 
	Compile in vim on linux mint
	:w | ! gcc qh_1.c -lm && ./a.out

	Sample run.
/usr/bin/ld: warning: /tmp/cctdQsZ5.o: requires executable stack (because the .note.GNU-stack section is executable)
1, f1, -1.110223024625157e-16
2, f2, 2.220446049250313e-16
3, f3, -6.661338147750939e-16
4, f4, 6.661338147750939e-16
5, f5, 6.661338147750939e-16
6, f6, 0.000000000000000e+00
7, f5, 6.661338147750939e-16
8, f5, -6.661338147750939e-16
9, f6, -1.665334536937735e-16
10, f6, 2.220446049250313e-16

*/
