function Qh(f,a=0,b=1,h=1/64,wanna_print=false)
  E(x)=exp.(-exp.(x.*(1+exp.(-x))))
  Ep(x)=exp.(-exp.(x.*(1+exp.(-x)))).* exp.(x.*(1+exp.(-x))).* (1+(1-x).*exp.(-x))
  if a==b
    return 0
  end
  if a>b
    return Qh(x->-f(x),b,a,h,wanna_print)
  end
  if a==-Inf && b==Inf
    return Qh(x->f(x/(1-x^2))*(1+x^2)/(1-x^2)^2,-1,1,h,wanna_print)
  end
  if isfinite(a) && b==Inf
    return Qh(x->f(a-1+1/x)/x^2,0,1,h,wanna_print)
  end
  if a==-Inf && isfinite(b)
    return Qh(x->f(b+1-1/x)/x^2,0,1,h,wanna_print)
  end
  ba=b-a
  id(x)=ifelse.(isfinite.(x),x,0)
  g(x)=ba*id(f(a+ba*E(x)).*Ep(x))
  S=h*g(0)
  n=-1 ; x=-h ; SS=S+h*g(x)
  while SS!=S || n>-10
    n -= 1 ; x -= h
    S=SS
    SS += h*g(x)
  end
  N=-n
  n=1 ; x=h ; SS=S+h*g(x)
  while SS!=S || n<10
    n += 1 ; x += h
    S=SS
    SS += h*g(x)
  end
  N += n
  if wanna_print
    println("N=",N)
  end
  return S
end

#=
Qh(x->1/(1+x^2),0,1,2^-3)-pi/4

Qh(x->exp(-pi*(x/(1-x^2))^2)*(1+x^2)/(1-x^2)^2,-1,1)

xi=rand();Qh(x->exp(-pi*x^2)*exp(2im*pi*xi*x),-Inf,Inf,2^-6)-exp(-pi*xi^2)

Qh(x->1/((x-1)^2+1),-Inf,Inf)-pi

Qh(x->Qh(y->1/(1-x*y),0,1),0,1)-pi^2/6

Qh(x->Qh(y->Qh(z->exp(-pi*(x^2+y^2+z^2)),-Inf,Inf),-Inf,Inf),-Inf,Inf)

Qh(x->1/(x^2+x^4+x^6+1),-Inf,Inf)-pi/2

Qh(x->Qh(y->Qh(z->1/(1+x^8+y^8+z^8),-Inf,Inf),-Inf,Inf),-Inf,Inf)

Qh(x->Qh(y->Qh(z->1/(1+x^8+y^8+z^8),-Inf,Inf,2^-7),-Inf,Inf,2^-7),-Inf,Inf,2^-7)

J=1.312142259171672454121313174812588624515937594782958395751511456250238626552831807801473842246748085 ;
for i=4:20; println(Qh(x->exp(sin(50*cos(60*sin(7*x)))),0,1,2.0^-i)-J);end

Qh(x->sinc(x)^2,0,Inf,2^-8)

Qh(x->sinc(x)^2,0,Inf,2^-10)

@time Qh(x->Qh(y->Qh(z->1/(1+x^2+y^2+z^2)^2,-Inf,Inf),-Inf,Inf),-Inf,Inf)-pi^2

	Results.
julia> Qh(x->1/(1+x^2),0,1,2^-3)-pi/4
-1.1102230246251565e-16
julia> Qh(x->exp(-pi*(x/(1-x^2))^2)*(1+x^2)/(1-x^2)^2,-1,1)
0.9999999999999999
julia> xi=rand();Qh(x->exp(-pi*x^2)*exp(2im*pi*xi*x),-Inf,Inf,2^-6)-exp(-pi*xi^2)
-6.938893903907228e-18 - 7.121435436161413e-17im
julia> Qh(x->1/((x-1)^2+1),-Inf,Inf)-pi
1.7763568394002505e-15
julia> Qh(x->Qh(y->1/(1-x*y),0,1),0,1)-pi^2/6
2.220446049250313e-16
julia> Qh(x->Qh(y->Qh(z->exp(-pi*(x^2+y^2+z^2)),-Inf,Inf),-Inf,Inf),-Inf,Inf)
1.0000000000000002
julia> Qh(x->1/(x^2+x^4+x^6+1),-Inf,Inf)-pi/2
2.220446049250313e-16
julia> Qh(x->Qh(y->Qh(z->1/(1+x^8+y^8+z^8),-Inf,Inf),-Inf,Inf),-Inf,Inf)
9.5850202080102
julia> Qh(x->Qh(y->Qh(z->1/(1+x^8+y^8+z^8),-Inf,Inf,2^-7),-Inf,Inf,2^-7),-Inf,Inf,2^-7)
9.58502020801018
julia> J=1.312142259171672454121313174812588624515937594782958395751511456250238626552831807801473842246748085 ;
julia> for i=4:20; println(Qh(x->exp(sin(50*cos(60*sin(7*x)))),0,1,2.0^-i)-J);end
-0.13861826382441578
-0.04491236394388798
-0.05691440171122775
-0.011401736434088816
-0.02039126438055061
-0.006064526490320743
-0.02113910156564991
-0.013437778777823484
-0.0006958604690183723
0.00011286802666243823
1.3414691579782811e-10
-1.2501111257279263e-13
-2.262634524186069e-13
-5.588862705963038e-13
-1.1377565556358604e-12
-2.2799540033702215e-12
-4.4790837705477315e-12
julia> Qh(x->sinc(x)^2,0,Inf,2^-8)
0.4999864674539415
julia> Qh(x->sinc(x)^2,0,Inf,2^-10)
0.4999930545694503
julia> @time Qh(x->Qh(y->Qh(z->1/(1+x^2+y^2+z^2)^2,-Inf,Inf),-Inf,Inf),-Inf,Inf)-pi^2
  7.322418 seconds (7.31 M allocations: 313.557 MiB, 1.19% gc time, 19.89% compilation time)
-8.881784197001252e-15
=#

