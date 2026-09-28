% u_exaxct: exact solution
%function [ result ] = u_exact(x, alpha,t1,beta)
function [ result ] = u_exact(x,alpha)
%mu=1;a1=2;b1=beta;
%b1=4/3;
% if ((sum(x.^2))>1)
%     result = 0;
% else
 %result = ml(-t1.^beta,beta,1).*(1-sum(x.^2)).^(alpha/2);
 %result=(t1.^a1).*exp(-1*sum(x.^2));
 %result = (t1.^(a1)).*(1-sum(x.^2)).^(alpha/2);
 %result = (1-sum(x.^2)).^(alpha/2);
 result=x.^2;
 %result=-1/4*(x(1).^2+x(2).^2);
 %result = (t1.^(a1)).*(1+sum(x.^2)).^(-7/2);
 %result = (t1.^(a1).*ml(1*t1.^mu,mu,a1+1)+t1.^b1).*(1-sum(x.^2)).^(alpha/2);
 %result = t1.^(a1).*(1-sum(x.^2)).^(alpha/2)*(2+alpha/2)/2*((alpha/2+3)*sum(x.^2)-2)/sqrt(pi).*(x(1,:));
end
