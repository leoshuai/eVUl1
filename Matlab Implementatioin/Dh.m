function g= Dh( x )
n=length(x);
g=zeros(n,1);
t=ones(n,1);
g(x>0)=t(x>0);
g(x<0)=-t(x<0);
end

