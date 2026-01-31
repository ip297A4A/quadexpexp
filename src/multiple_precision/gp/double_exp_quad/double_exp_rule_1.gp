Qh(f,a=0,b=1,h=1/64.0,wanna_print=1,K=10)={
  my(E,Ep,ba,g,S,N,n,x);
  (E=x->exp(-exp(x*(1+exp(-x))))) ;
  (Ep=x->exp(-exp(x*(1+exp(-x))))*exp(x*(1+exp(-x)))*(1+(1-x)*exp(-x))) ;
  if(a==b,
    return(0)) ;
  if(a>b,
    return(-Qh(f,b,a,h,wanna_print,K))) ;
  if(a==-oo && b==+oo,
    return(Qh((x->f(x/(1-x^2))*(1+x^2)/(1-x^2)^2),-1,1,h,wanna_print,K))) ;
  if((-oo<a)&&(a<+oo) && b==+oo,
    return(Qh((x->f(a-1+1/x)/x^2),0,1,h,wanna_print,K))) ;
  if(a==-oo && (-oo<b)&&(b<+oo),
    return(Qh((x->f(b+1-1/x)/x^2),0,1,h,wanna_print,K))) ;
  ba=b-a ;
  (id=x->if((-oo<x)&&(x<+oo),x,0)) ;
  (g=x->ba*id(f(a+ba*E(x))*Ep(x))) ;
  S=h*g(0) ;
  n=-1 ; x=- h ; SS=S+h*g(x) ;
  while(SS!=S || n>-K,
    n=n-1 ; x=x-h ; S=SS ; SS=SS+h*g(x) ;
  ) ;
  N=-n ;
  n=1 ; x=h ; SS=S+h*g(x) ;
  while(SS!=S || n<K,
    n=n+1 ; x=x+h ; S=SS ; SS=SS+h*g(x) ;
  ) ;
  N=N+n+1 ;
  if(wanna_print, print("N=",N,"\n"));
  return(S) ;
}

/*	Sample run.
\p100
#
Qh(x->1/(1+x^2),-oo,oo,2^-8)-Pi
Qh(x->1/sqrt(x),0,1,2^-8)-2
Qh(x->sqrt(1-x^2),-1,1,2^-8)-Pi/2
Qh(x->1/sqrt(1-x^2),-1,1,2^-8)-Pi
Qh(x->-log(-log(x)),0,1,2^-8)-Euler
Qh(x->exp(-x)*log(x),0,oo,2^-8)+Euler
2*Qh(x->(exp(-x^2)-exp(-x))/x,0,oo,2^-9)-Euler
Qh(x->exp(-Pi*x^2),-oo,oo,2^-9)-1
xi=(random(1.0)-0.5)+I*(random(1.0)-0.5);abs(Qh(x->exp(-Pi*x^2)*exp(-2*Pi*I*xi*x),-oo,oo,2^-9)/exp(-Pi*xi^2)-1)
Z(s,h)=1/(s-1)+1/2+2*Qh(t->sin(s*atan(t))/(1+t^2)^(s/2)/(exp(2*Pi*t)-1),0,oo,h)
Z(2,2^-8)-Pi^2/6
Z(4,2^-8)-Pi^4/90
Z(-1,2^-8)-(-1/12)
Z(-3,2^-8)-1/120
s=(10*(-0.5+random(1.0)))+I*(10*(-0.5+random(1.0)));abs(Z(s,2^-8)/zeta(s)-1)
t=1e-50;(Z(1+t,2^-8)+Z(1-t,2^-8))/2-Euler
Qh(x->exp(-x)*cos(x),0,oo,2^-9)-1/2
Qh(x->exp(-x)*sin(x),0,oo,2^-9)-1/2
Qh(exp,-oo,1,2^-8)-exp(1)

	Results.
(14:22) gp >
(14:23) gp > \p100
   realprecision = 115 significant digits (100 digits displayed)
(14:23) gp > #
   timer = 1 (on)
(14:23) gp > Qh(x->1/(1+x^2),-oo,oo,2^-8)-Pi
N=2489

time = 187 ms.
%3 = -9.694937818545780033 E-113
(14:23) gp > Qh(x->1/sqrt(x),0,1,2^-8)-2
N=2667

time = 188 ms.
%4 = -5.121566627702980133 E-113
(14:23) gp > Qh(x->sqrt(1-x^2),-1,1,2^-8)-Pi/2
N=2298

time = 156 ms.
%5 = -5.405816113482332846 E-113
(14:23) gp > Qh(x->1/sqrt(1-x^2),-1,1,2^-8)-Pi
N=2813

time = 203 ms.
%6 = -1.0953756969854342048 E-112
(14:23) gp > Qh(x->-log(-log(x)),0,1,2^-8)-Euler
N=2493

time = 220 ms.
%7 = 5.913404480945462690 E-114
(14:23) gp > Qh(x->exp(-x)*log(x),0,oo,2^-8)+Euler
N=1420

time = 125 ms.
%8 = 7.055478307737504841 E-114
(14:23) gp > 2*Qh(x->(exp(-x^2)-exp(-x))/x,0,oo,2^-9)-Euler
N=2827

time = 219 ms.
%9 = -1.4187094870594479158 E-113
(14:23) gp > Qh(x->exp(-Pi*x^2),-oo,oo,2^-9)-1
N=978

time = 78 ms.
%10 = -1.1166944084188856582 E-113
<dom(1.0)-0.5);abs(Qh(x->exp(-Pi*x^2)*exp(-2*Pi*I*xi*x),-oo,oo,2^-9)/exp(-Pi*xi^2)-1)
N=979

time = 94 ms.
%11 = 2.3146265080628421821 E-114
(14:23) gp > Z(s,h)=1/(s-1)+1/2+2*Qh(t->sin(s*atan(t))/(1+t^2)^(s/2)/(exp(2*Pi*t)-1),0,oo,h)
%12 = (s,h)->1/(s-1)+1/2+2*Qh(t->sin(s*atan(t))/(1+t^2)^(s/2)/(exp(2*Pi*t)-1),0,oo,h)
(14:23) gp > Z(2,2^-8)-Pi^2/6
N=1308

time = 109 ms.
%13 = -3.502359735495595928 E-114
(14:23) gp > Z(4,2^-8)-Pi^4/90
N=1306

time = 110 ms.
%14 = -2.1318711433451453476 E-114
(14:23) gp > Z(-1,2^-8)-(-1/12)
N=1312

time = 126 ms.
%15 = 1.7786742376520508307 E-114
(14:23) gp > Z(-3,2^-8)-1/120
N=1313

time = 109 ms.
%16 = 2.1855909122349932562 E-114
(14:23) gp > s=(10*(-0.5+random(1.0)))+I*(10*(-0.5+random(1.0)));abs(Z(s,2^-8)/zeta(s)-1)
N=1313

time = 156 ms.
%17 = 8.905031169727938548 E-114
(14:23) gp > t=1e-50;(Z(1+t,2^-8)+Z(1-t,2^-8))/2-Euler
N=1310

N=1310

time = 266 ms.
%18 = -1.1665312939063460046 E-85
(14:23) gp > Qh(x->exp(-x)*cos(x),0,oo,2^-9)-1/2
N=2828

time = 203 ms.
%19 = -2.0430431790390976247 E-113
(14:23) gp > Qh(x->exp(-x)*sin(x),0,oo,2^-9)-1/2
N=2541

time = 203 ms.
%20 = -1.0735493971845196214 E-113
(14:23) gp > Qh(exp,-oo,1,2^-8)-exp(1)
N=1415

time = 110 ms.
%21 = -6.497131103528062012 E-113
(14:23) gp > \p1001
   realprecision = 1001 significant digits
(14:23) gp > round(-log(abs(  Qh(x->1/(1+x^2),0,1,2^-12)/(Pi/4)-1  ))/log(10))
N=56103

time = 1min, 39,923 ms.
%22 = 997
(14:25) gp >
*/
