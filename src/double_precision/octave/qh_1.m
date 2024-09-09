function [S,N]=Qh(f,a,b,h=1/64)
  E=@(x)exp(-exp(x.*(1+exp(-x)))) ;
  Ep=@(x)exp(-exp(x.*(1+exp(-x)))).* exp(x.*(1+exp(-x))).* (1+(1-x).*exp(-x)) ;
  if a==b
    S=0 ; N=0 ; return ;
  end
  if a==-inf && b==inf
    [S,N]=Qh(@(x)f(x./(1-x.^2)).*(1+x.^2)./(1-x.^2).^2,-1,1,h) ;
    return ;
  end
  if isfinite(a) && b==inf
    [S,N]=Qh(@(x)f(a+x./(1-x))./(1-x).^2,0,1,h) ;
    return ;
  end
  if a==-inf && isfinite(b)
    [S,N]=Qh(@(x)f(b-x./(1-x))./(1-x).^2,0,1,h) ;
    return ;
  end
  ba=b-a ;
  id=@(x)merge(isfinite(x),x,0) ;
  g=@(x)ba*id(f(a+ba*E(x)).*Ep(x)) ;
  S=h*g(0) ;
  n=-1 ; x=-h ; SS=S+h*g(x) ;
  while SS!=S || n>-10
    n-- ; x -= h ;
    S=SS ;
    SS += h*g(x) ;
  end
  N=-n ;
  n=1 ; x=h ; SS=S+h*g(x) ;
  while SS!=S || n<10
    n++ ; x += h ;
    S=SS ;
    SS += h*g(x) ;
  end
  N += n+1 ;
end


function [S,N]=Qts(f,a,b,h=1/64)
  E=@(x)tanh(pi/2*sinh(x)) ;
  Ep=@(x)sech(pi/2*sinh(x)).^2.*cosh(x).*pi/2 ;
  if a==b
    S=0 ; N=0 ; return ;
  end
  if a==-inf && b==inf
    [S,N]=Qts(@(x)f(x./(1-x.^2))*(1+x.^2)/(1-x.^2).^2,-1,1,h) ;
    return ;
  end
  if isfinite(a) && b==Inf
    [S,N]=Qts(@(x)f(a+x./(1-x))./(1-x).^2,0,1,h) ;
    return ;
  end
  if a==-Inf && isfinite(b)
    [S,N]=Qts(@(x)f(b-x./(1-x))./(1-x).^2,0,1,h) ;
    return ;
  end
  ba=(b-a)/2 ; m=(a+b)/2 ;
  id=@(x)merge(isfinite(x),x,0) ;
  g=@(x)ba*id(f(m+ba*E(x)).*Ep(x)) ;
  S=h*g(0) ;
  n=-1 ; x=-h ; SS=S+h*g(x) ;
  while SS!=S || n>-10
    n-- ; x -= h ;
    S=SS ;
    SS += h*g(x) ;
  end
  N=-n ;
  n=1 ; x=h ; SS=S+h*g(x) ;
  while SS!=S || n<10
    n++ ; x += h ;
    S=SS ;
    SS += h*g(x) ;
  end
  N += n+1 ;
end

%{
>> Qh(f=@(x)exp(sin(11*x.^2))+cos(19*exp(sin(13*x))),a=0,b=1,2^-12)-1.4999465049014194206968363892452710985
ans = -1.154631945610163e-14
>> Qts(f=@(x)exp(sin(11*x.^2))+cos(19*exp(sin(13*x))),a=0,b=1,2^-12)-1.4999465049014194206968363892452710985
ans = -6.883382752675971e-15
>> Qh(f=@(x)sinc(x),a=0,b=1000,2^-9)
ans = 0.499898678836888
>> Qh(f=@(x)sinc(x),a=0,b=1000,2^-9)-0.49989867883688960177412531246608695583
ans = -2.053912595556540e-15
>> Qts(f=@(x)sinc(x),a=0,b=1000,2^-9)-0.49989867883688960177412531246608695583
ans = -4.551914400963142e-15
>> Qh(@(x)log(2+cos(x+x.^2)),0,100,2^-15)-62.518029770285939532930943231054768985
ans = -6.899369964230573e-12
>> Qts(@(x)log(2+cos(x+x.^2)),0,100,2^-15)-62.518029770285939532930943231054768985
ans = -3.048228336410830e-12
>> Qh(@exp,-inf,1)-exp(1)
ans = -1.776356839400250e-15
>> Qts(@exp,-inf,1)-exp(1)
ans = -1.332267629550188e-15
>> Qh(@(x)exp(-pi*x.^2),-inf,inf)-1
ans = -1.110223024625157e-16
>> Qts(@(x)exp(-pi*x.^2),-inf,inf)-1
ans = -1.110223024625157e-16
>> Qh(@(x)1./sqrt(x),0,1)-2
ans = 4.440892098500626e-16
>> Qts(@(x)1./sqrt(x),0,1)-2
ans = -1.061142551606054e-08
>> Qh(@(x)log(x),0,1)+1
ans = 4.440892098500626e-16
>> Qts(@(x)log(x),0,1)+1
ans = 2.109423746787797e-15
>>
%}
