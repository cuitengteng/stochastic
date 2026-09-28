%This script is for simulating the jump distance
function result= jumpdis(rand_omega,alpha,r)


% rand_omega: rand parameter
% r the radius of small ball
% alpha: fractional index

k = pi*rand_omega / ( sin(pi*alpha/2)*beta(1-alpha/2,alpha/2) );
result = r / sqrt(1-betaincinv(k,1-(alpha/2),alpha/2));
end
