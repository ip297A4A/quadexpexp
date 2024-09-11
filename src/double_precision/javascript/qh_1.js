Qh=function(f,a,b,h=1/64){
  var exp=Math.exp ;
  var E=function(x){return(exp(-exp(x*(1+exp(-x))))) ;} ;
  var Ep=function(x){return(exp(-exp(x*(1+exp(-x))))* exp(x*(1+exp(-x)))* (1+(1-x)*exp(-x))) ;} ;
  if(a==b){return(0) ; } 
  if(a>b){return(-Qh(f,b,a,h)); }
  if(a==-Infinity && b==Infinity){
    return(Qh(x=> f(x/(1-x**2))*(1+x**2)/(1-x**2)**2,-1,1,h)) ;
  }
  if(isFinite(a) && b==Infinity){
    return(Qh(x=> f(a-1+1/x)/x**2,0,1,h)) ;
  }
  if(a==-Infinity && isFinite(b)){
    return(Qh(x=> f(b+1-1/x)/x**2,0,1,h)) ;
  }
  var ba=b-a ;
  var id=function(x){return(isFinite(x) ? x : 0) ;} ;
  var g=function(x){return(ba*id(f(a+ba*E(x))*Ep(x))) ;} ;
  var S=h*g(0) ;
  var n=-1 ; var x=-h ; var SS=S+h*g(x) ;
  while(SS!=S || n>-10){
    n-- ; x -= h ; S=SS ; SS += h*g(x) ;
  }
  var N=-n ; n=1 ; x=h ; SS=S+h*g(x) ;
  while(SS!=S || n<10){
    n++ ; x += h ; S=SS ; SS += h*g(x) ;
  }
  N += n+1 ;
  console.log(N) ;
  return(S) ;
}

