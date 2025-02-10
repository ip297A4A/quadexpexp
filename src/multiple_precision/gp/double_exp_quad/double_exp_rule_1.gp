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
    return(Qh((x->f(a+x/(1-x))/(1-x)^2),0,1,h,wanna_print,K))) ;
  if(a==-oo && (-oo<b)&&(b<+oo),
    return(Qh((x->f(b-x/(1-x))/(1-x)^2),0,1,h,wanna_print,K))) ;
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
(18:20) gp > \p100
(18:20) gp > #
   timer = 1 (on)
(18:20) gp > Qh(x->1/(1+x^2),-oo,oo,2^-8)-Pi
N=2489

time = 187 ms.
%302 = -9.694937818545780033 E-113
(18:20) gp > Qh(x->1/sqrt(x),0,1,2^-8)-2
N=2667

time = 188 ms.
%303 = -5.121566627702980133 E-113
(18:20) gp > Qh(x->sqrt(1-x^2),-1,1,2^-8)-Pi/2
N=2298

time = 172 ms.
%304 = -5.405816113482332846 E-113
(18:20) gp > Qh(x->1/sqrt(1-x^2),-1,1,2^-8)-Pi
N=2813

time = 203 ms.
%305 = -1.0953756969854342048 E-112
(18:20) gp > Qh(x->-log(-log(x)),0,1,2^-8)-Euler
N=2493

time = 251 ms.
%306 = 5.913404480945462690 E-114
(18:20) gp > Qh(x->exp(-x)*log(x),0,oo,2^-8)+Euler
N=1751

time = 172 ms.
%307 = 1.4973856840162330417 E-113
(18:20) gp > 2*Qh(x->(exp(-x^2)-exp(-x))/x,0,oo,2^-9)-Euler
N=3489

time = 312 ms.
%308 = -3.253641435438662304 E-113
(18:20) gp > Qh(x->exp(-Pi*x^2),-oo,oo,2^-9)-1
N=978

time = 78 ms.
%309 = -1.1166944084188856582 E-113
<dom(1.0)-0.5);abs(Qh(x->exp(-Pi*x^2)*exp(-2*Pi*I*xi*x),-oo,oo,2^-9)/exp(-Pi*xi^2)-1)
N=982

time = 94 ms.
%310 = 7.742105569865964861 E-114
(18:20) gp > Z(s,h)=1/(s-1)+1/2+2*Qh(t->sin(s*atan(t))/(1+t^2)^(s/2)/(exp(2*Pi*t)-1),0,oo,h)
%311 = (s,h)->1/(s-1)+1/2+2*Qh(t->sin(s*atan(t))/(1+t^2)^(s/2)/(exp(2*Pi*t)-1),0,oo,h)
(18:20) gp > Z(2,2^-8)-Pi^2/6
N=1684

time = 189 ms.
%312 = -4.009948102958725773 E-114
(18:20) gp > Z(4,2^-8)-Pi^4/90
N=1683

time = 172 ms.
%313 = -4.212983449943977711 E-114
(18:20) gp > Z(-1,2^-8)-(-1/12)
N=1686

time = 172 ms.
%314 = 2.1107216280341816040 E-114
(18:20) gp > Z(-3,2^-8)-1/120
N=1687

time = 171 ms.
%315 = 4.209599527494223512 E-114
(18:20) gp > s=(10*(-0.5+random(1.0)))+I*(10*(-0.5+random(1.0)));abs(Z(s,2^-8)/zeta(s)-1)
N=1685

time = 219 ms.
%316 = 8.629449085257190386 E-114
(18:20) gp > t=1e-50;(Z(1+t,2^-8)+Z(1-t,2^-8))/2-Euler
N=1685

N=1685

time = 361 ms.
%317 = -1.1665312939063460046 E-85
(18:20) gp > Qh(x->exp(-x)*cos(x),0,oo,2^-9)-1/2
N=3489

time = 312 ms.
%318 = -1.7765592861209544563 E-113
(18:20) gp > Qh(x->exp(-x)*sin(x),0,oo,2^-9)-1/2
N=3124

time = 281 ms.
%319 = -1.6369724850685937490 E-113
(18:20) gp > Qh(exp,-oo,1,2^-8)-exp(1)
N=1745

time = 157 ms.
%320 = -7.898074997726300383 E-113
(18:20) gp >


	Another run.
(17:08) gp > \p1001
(17:08) gp > round(-log(abs(  Qh(x->1/(1+x^2),0,1,2^-12)/(Pi/4)-1  ))/log(10))
N=56103

time = 1min, 44,361 ms.
%35 = 997
(17:11) gp >
*/
