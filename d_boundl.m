% function d = d_boundl(x, g)
%     % 判断点是否在局部区域内部（三维）
%     % x: 3x1 向量 [x; y; z]
%     % g: 曲面函数 x = g(y, z)
%     % 返回: 1 表示在局部区域内部，0 表示在外部
% 
%     xc = x(1);
%     yc = x(2);
%     zc = x(3);
% 
%     % 曲面上的x坐标
%     x_interface = g(yc, zc);
% 
%     % 局部区域：x < g(y, z)
%     is_in_local = (xc <= x_interface);
% 
%     % 同时检查是否在立方体边界内
%     is_in_cube = (xc >= -1 && xc <= 1 && yc >= -1 && yc <= 1 && zc >= -1 && zc <= 1);
% 
%     d = is_in_local && is_in_cube;
% end

function [ d ] = d_boundl( x )
d = 0;  % 默认返回值

    d_left = x(1) - (-1);      % 到左边界 x = -1 的距离
    d_right = 0 - x(1);        % 到右边界 x = 0 的距离
    d_bottom = x(2) - (-1);    % 到下边界 y = -1 的距离
    d_top = 1 - x(2);          % 到上边界 y = 1 的距离

    % 取最小值
    d = max(min([d_left, d_right, d_bottom, d_top]),0); 



