%%%% Subroutine for computing the exit time in time
function result = rand_stable(alpha, theta, t)
    % 生成一个随机数
  
    random_value = rand();
    expon_value = exprnd(1); % MATLAB 中的 exponential random variable 生成
    % 计算 sigma 函数
    if random_value ==0
    result= (theta * t).^(1 / alpha).* (1 - alpha).^((1 - alpha) / alpha).*alpha ./ (expon_value).^((1 - alpha) / alpha);
      else 
     result=(theta * t).^(1 / alpha).* (sin((1 - alpha) * pi * random_value)).^((1 - alpha) / alpha) ...
                               .* (sin(alpha * pi * random_value)) ...
                               ./ (sin(pi * random_value)).^(1 / alpha)./ (expon_value).^((1 - alpha) / alpha);
    end
    %sigma_value = sigma(alpha, random_value);
    %sigma_value=1;
    % 生成指数分布随机数
