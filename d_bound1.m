function [ d ] = d_bound1( x )
d=0;
if x(1)>-1&&x(1)<=0&&x(2)>-1&&x(2)<1
    d_left = x(1) - (-1);      % 到左边界 x = -1 的距离
    d_right = 0 - x(1);        % 到右边界 x = 0 的距离
    d_bottom = x(2) - (-1);    % 到下边界 y = -1 的距离
    d_top = 1 - x(2);          % 到上边界 y = 1 的距离
  d = max(min([d_left, d_right, d_bottom, d_top]),0);  
    % 取最小值
elseif x(1)>0&&x(1)<1&&x(2)>-1&&x(2)<1
    d_left = x(1) - (0);      % 到左边界 x = -1 的距离
    d_right = 1- x(1);        % 到右边界 x = 0 的距离
    d_bottom = x(2) - (-1);    % 到下边界 y = -1 的距离
    d_top = 1 - x(2);          % 到上边界 y = 1 的距离
    
    % 取最小值
    d = max(min([d_left, d_right, d_bottom, d_top]),0);

end