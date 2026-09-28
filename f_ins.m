
% f_ins: source term 

%function [ result ] = f_ins(x,t,alpha,n,beta)
function [ result ] = f_ins(x,alpha)
%mu=1;a1=2;b1=beta;
%c1 = [pi/3, -pi/4];
%c2 = [-pi/2, 2*pi/3];
a1=1/3;
c1=pi/6;
c2=pi/3;
n=2;
%b1=4/3;
%result=zeros(10,11);
% if (sqrt(sum(x.^2))>=1)
%      result = 0;
% else
%result=-ml(-t.^beta,beta,1).*exp(-1*sum(x.^2))+ml(-t.^beta,beta,1).*kummerM_series((n+alpha)/2,n/2,-1*sum(x.^2),500,1e-16);
%result =t.^(a1).* 2^(alpha)*gamma(n/2+alpha/2)/gamma(n/2).*hypergeom((n+alpha)/2,n/2,-1*sum(x.^2))+exp(-1*sum(x.^2))*gamma(a1+1)/gamma(a1+1-beta)*t.^(a1-beta);
 %result =(1-sum(x.^2)).^(alpha/2)*(-sin(t)./(1+10*t.^2)-20*t.*cos(t)./(1+10*t.^2).^2 )+cos(t)./(1+10*t.^2)*(2^alpha)*gamma(1+(alpha/2))*gamma(n/2+(alpha/2))/gamma(n/2);
 %result =(1-sum(x.^2)).^(alpha/2)*(gamma(a1+1)/gamma(a1+1-beta)*t.^(a1-beta))+(t.^(a1)).*(2^alpha)*gamma(1+(alpha/2))*gamma(n/2+(alpha/2))/gamma(n/2);
 %result =(2^alpha)*gamma(1+(alpha/2))*gamma(n/2+(alpha/2))/gamma(n/2);
 %result=t^(1/3).*cos(x(1)*x(2));
 result=1+1/2*cos(x(1)*x(2));
 %result=1+1/2*x(1).*x(2);
 %result =(1+sum(x.^2)).^(-7/2)*(gamma(a1+1)/gamma(a1+1-beta)*t.^(a1-beta))+(t.^(a1))*2^alpha*gamma(alpha/2+7/2)*gamma(alpha/2+1)/gamma(7/2)*hypergeom_transformed((2+alpha)/2,(7+alpha)/2,n/2, -sum(x.^2),1e-16);
 %result =(1-sum(x.^2)).^(alpha/2).*(t.^(a1-beta).*ml(1*t.^mu,mu,a1+1-beta)+gamma(b1+1)/gamma(b1+1-beta)*t.^(b1-beta))+(t.^a1.*ml(1*t.^mu,mu,a1+1)+t.^b1).*(2^alpha)*gamma(1+(alpha/2))*gamma(n/2+(alpha/2))/gamma(n/2);
 %result = (gamma(a1+1)/gamma(a1+1-beta)*t.^(a1-beta)).*(1-sum(x.^2)).^(alpha/2)*(2+alpha/2)/2*((alpha/2+3)*sum(x.^2)-2)/sqrt(pi)*(x(1,:))+...
     %t.^(a1).*(2^alpha)/4*(gamma(3+alpha/2))^2*((alpha/2+3)*sum(x.^2)-2)/sqrt(pi).*(x(1,:));
%result= t.^(a1).*(cos(c2*x(1,:)).^(2) + sin(c1*x(2,:)).^(2)-(alpha*x(1,:)*x(2,:))^3);
%result= (cos(c2*x(1,:)^2-x(1,:)*x(2,:)) + sin(c1*x(2,:)^2+x(1,:)*x(2,:)));
end

