%This script is for simulating the jump distance
function result= H(y,alpha,r)
% if alpha < 3e-1
%         result = r * sqrt(alpha / (2 * y));
%         return;
% end
% k = pi*y/(sin(pi*alpha/2)*beta(1-alpha/2,alpha/2));
% result = r/(sqrt(1-betaincinv(k,1-(alpha/2),alpha/2)));   
result = r/(sqrt(real(betaincinv(y,(alpha/2),1-alpha/2))));
end
