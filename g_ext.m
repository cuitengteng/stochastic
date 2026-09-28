
% g_ext: exterior function

function [ result ] = g_ext( x,t1,alpha,beta )
a1=2;
%if (sum(x.^2))<=1
result=0;
%result=t1.^a1.*exp(-1*sum(x.^2));
%result = (t1.^(a1)).*(1+sum(x.^2)).^(-7/2);
%result = (t1.^(a1)).*(1-sum(x.^2)).^(alpha/2);
%result =exp(-1*sum(x.^2))*t1.^(a1);
end

