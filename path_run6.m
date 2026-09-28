%function [error] = path_run6(alpha,T,dt,num_path)%%%%%main Steady-state equation
% Omega = [-1, 1];                % 全局域
% Omega_l = [-1, 0];          % 局部区域 (布朗运动)
% Omega_nl = [0 1]; % 非局部区域 (跳跃)
lx=-1;hx=1;ly=-1;hy=1;
num_mesh_x = 20; % 
refN=100;
num_mesh_y=num_mesh_x;
h = 2.0 / (num_mesh_x - 1);
x = linspace(lx, hx, num_mesh_x);
y = linspace(ly, hy, num_mesh_y);

[X, Y] = meshgrid(x, y);
x0=[(X(:))';(Y(:))'];
alpha1=2;
N=num_mesh_x*num_mesh_y;
num_path =3200;
alpha =2;
dt= 1e-4;
epoch=floor(num_path/dt);
%num=0;

n=2;
r=(dt*(2^alpha)*(gamma(1+(alpha/2)))^2)^(1/alpha);
r1=(dt*(2^alpha1)*(gamma(1+(alpha1/2)))^2)^(1/alpha1);
c1=(2^(alpha-1)*alpha*gamma(1+alpha/2)/pi/gamma(1-alpha/2));
%T =100;

solution = zeros(N,1);
tic;
if alpha == 2
    % alpha=2：左侧扩散系数 1，右侧 2，界面满足连续性与通量连续。
    M_list = [100 200 400 800 1600 3200];
    M_list = unique([M_list(M_list <= num_path), num_path]);
    [solution, solution_M] = wos_alpha2(x0, dt, M_list, 1e-8, 20260909);
else
parfor k=1:N
    sol_stoc_new = 0; 
    num=0;
 for p=1: num_path 
     num=num+1;
     t0 = 0;
    x = x0(:,k);
  estimate=0;
 for i=2: epoch
    t0=t0+dt;
    if d_bound(x) <= 0
        break;
    elseif d_boundl(x)>0 || (x(1)==0&&x(2)>-1&&x(2)<1)
        % Corrected local-region path:
        % first Brownian endpoint step, then one auxiliary stable candidate.
        estimate= estimate+f_ins(x,alpha1)*dt;
        r_new = H(rand,alpha1,r1);
        theta_new = rand*2*pi;
        z = x + r_new.*[cos(theta_new); sin(theta_new)];

        % If Brownian endpoint leaves Omega, terminate this path.
        if d_bound(z) <= 0
            break;


        % If Brownian endpoint crosses the flat interface x=0, reflect it
        % back to Omega_l.  Here this is z_1 -> -|z_1|.
        elseif z(1)>=0 && z(1)<1 && z(2)>-1 && z(2)<1
            z(1) = -abs(z(1));
        end
        

        % From the accepted Brownian endpoint z, independently generate one
        % stable candidate and accept only Omega_l -> Omega_nl jumps.
        r_new = H(rand,alpha,r);
        theta_new = rand*2*pi;
        y_cand = z + r_new.*[cos(theta_new); sin(theta_new)];

        if d_boundnl(y_cand)>0
            x = y_cand;
        else
            x = z;
        end
        


    elseif d_boundnl(x)>0
        % Corrected nonlocal-region path:
        % First complete jump B (Brownian-radius jump when alpha=2).
        % If B is still in Omega_nl, generate
        % a second independent candidate and accept it only inside Omega_nl.
        % Two accepted interior Brownian steps give twice the covariance.
        % Add the source only once per dt, including when both jumps occur.
        estimate=estimate+f_ins(x,alpha)*dt;
        r_new = H(rand,alpha,r);
        theta_new = rand*2*pi;
        b = x + r_new.*[cos(theta_new); sin(theta_new)];

        if d_bound(b) <= 0
            break;
        elseif d_boundl(b)>0
            x = b;
        elseif d_boundnl(b)>0

            r_new = H(rand,alpha,r);
            theta_new = rand*2*pi;
            y_cand = b + r_new.*[cos(theta_new); sin(theta_new)];

            if d_boundnl(y_cand)>0
                x = y_cand;
            else
                x = b;
            end
        else
            break;
        end

    else
        break;
    end
 end
if isnan(estimate) || isinf(estimate)
     estimate=0;
     num=num-1;
end              
sol_stoc_new=sol_stoc_new + estimate;

end
solution_stoc=sol_stoc_new/num;
solution(k)= solution_stoc;
end
end
u_exact=solve_coupling_fixed(alpha,refN,4096);
%u_exact=coupled_fractional_solver(alpha,refN);
solution_exact=u_exact;
%error1=abs(sqrt(1/N*sum((u_num_2d(:) -u_ref_subsampled(:)).^2))/sqrt(1/N*sum((u_ref_subsampled(:)).^2)))
error=calc_error_and_plot(solution, solution_exact, num_mesh_x, refN)
toc;
% figure('Position', [100, 100, 800, 800]);
%      u=reshape(solution, num_mesh_x,num_mesh_y);
%      u_exact=reshape(solution_exact,100,100);
%     % --- Plot 1: 2D Color Map ---
%     subplot(2, 2, 1);
%     imagesc([-1, 1], [-1, 1], u); 
% 
%     set(gca, 'YDir', 'normal'); % 对应 origin='lower'
%     colormap('jet');
%     colorbar;
%     hold on;
%     xline(0, '--w', 'Interface', 'LineWidth', 2);
%     title(sprintf('Solution u(x,y) [s=%.1f]', alpha/2));
%     xlabel('x (i index)');
%     ylabel('y (j index)');
%     hold off;
%     subplot(2, 2, 2);
%     imagesc([-1, 1], [-1, 1], u_exact); 
% 
%     set(gca, 'YDir', 'normal'); % 对应 origin='lower'
%     colormap('jet');
%     colorbar;
%     hold on;
%     xline(0, '--w', 'Interface', 'LineWidth', 2);
%     title(sprintf('Solution uexact(x,y) [s=%.1f]', alpha/2));
%     xlabel('x (i index)');
%     ylabel('y (j index)');
%     hold off;
%     % --- Plot 2: Cross Section ---
%     subplot(2, 2, 3);
% 
%     % 在 Python 中取的是 u[:, mid_y] (固定 j, 变 i)
%     mid_y = floor(num_mesh_x / 2) + 1;
% 
%     % u 已经是 (i, j) 格式了，所以取第 mid_y 列即可
%     % 注意：Python的 mid_y 对应这里的列索引
%     u_cross = u(mid_y,:);
% 
%     x_axis = linspace(-1, 1, num_mesh_x);
% 
%     plot(x_axis, u_cross, 'b-o', 'MarkerSize', 4, 'DisplayName', 'Cross-section y=0');
%     hold on;
%     xline(0, '--r', 'DisplayName', 'Interface (x=0)');
% 
%     u_max = max(u(:));
%     text(-0.5, u_max*0.2, sprintf('Local Region\n(-Delta + Int)'), 'HorizontalAlignment', 'center');
%     text(0.5, u_max*0.2, sprintf('Nonlocal Region\n(Double Interaction)'), 'HorizontalAlignment', 'center');
% 
%     title('u cross Profile');
%     grid on;
%     set(gca, 'GridAlpha', 0.3);
%     legend('show');
%     hold off;
%     subplot(2, 2, 4);
% 
%     % 在 Python 中取的是 u[:, mid_y] (固定 j, 变 i)
%     mid_y = floor(100 / 2) + 1;
% 
%     % u 已经是 (i, j) 格式了，所以取第 mid_y 列即可
%     % 注意：Python的 mid_y 对应这里的列索引
%     u_cross = u_exact(mid_y,:);
% 
%     x_axis = linspace(-1, 1, 100);
% 
%     plot(x_axis, u_cross, 'b-o', 'MarkerSize', 4, 'DisplayName', 'Cross-section y=0');
%     hold on;
%     xline(0, '--r', 'DisplayName', 'Interface (x=0)');
% 
%     u_max = max(u_exact(:));
%     text(-0.5, u_max*0.2, sprintf('Local Region\n(-Delta + Int)'), 'HorizontalAlignment', 'center');
%     text(0.5, u_max*0.2, sprintf('Nonlocal Region\n(Double Interaction)'), 'HorizontalAlignment', 'center');
% 
%     title('uexact Cross-section Profile');
%     grid on;
%     set(gca, 'GridAlpha', 0.3);
%     legend('show');
%     hold off;

function [solution, solution_M] = wos_alpha2(x0, dt, M_list, tol, seed)
%WOS_ALPHA2 平面界面传输问题的二维 WoS 模拟。
% 左侧扩散系数为 1，右侧为 2，外边界值为零。
% dt 是单次出球的平均时间上限；tol 控制界面投影和边界终止。

num_path = max(M_list);
num_point = size(x0, 2);
solution = zeros(num_point, 1);
solution_M = zeros(num_point, numel(M_list));

parfor k = 1:num_point
    stream = RandStream('Threefry', 'Seed', seed);
    stream.Substream = k;

    x = repmat(x0(1,k), num_path, 1);
    y = repmat(x0(2,k), num_path, 1);
    estimate = zeros(num_path, 1);
    active = find(max(abs([x, y]), [], 2) < 1-tol);
    step = 0;

    while ~isempty(active)
        step = step + 1;
        if step > 1000000
            error('路径步数超过保护上限，不能将截断路径作为完整样本。');
        end

        x_old = x(active);
        y_old = y(active);
        on_interface = abs(x_old) <= tol;
        x_old(on_interface) = 0;
        num_active = numel(active);

        % 单侧球使用本侧扩散系数；界面球使用两侧平均系数。
        a = 1 + double(x_old > 0);
        a(on_interface) = 1.5;
        radius = min(sqrt(4*a*dt), min(1-abs(x_old), 1-abs(y_old)));
        radius(~on_interface) = min(radius(~on_interface), ...
                                    abs(x_old(~on_interface)));

        % 根据球内 Green 函数采样源项，权重为平均退出时间。
        r_source = radius .* sqrt(rand(stream,num_active,1) .* ...
                                  rand(stream,num_active,1));
        theta_source = 2*pi*rand(stream,num_active,1);
        x_source = x_old + r_source.*cos(theta_source);
        y_source = y_old + r_source.*sin(theta_source);
        source_value = f_ins([x_source.'; y_source.'], 2);
        estimate(active) = estimate(active) + ...
            radius.^2./(4*a) .* source_value(:);

        % 单侧球：均匀圆周出球。
        theta = 2*pi*rand(stream,num_active,1);
        dx = radius.*cos(theta);
        dy = radius.*sin(theta);

        % 界面中心球：右半圆概率 2/3，左半圆概率 1/3。
        num_interface = sum(on_interface);
        if num_interface > 0
            theta = pi*(rand(stream,num_interface,1)-0.5);
            side = 2*double(rand(stream,num_interface,1) < 2/3)-1;
            dx(on_interface) = side.*radius(on_interface).*cos(theta);
            dy(on_interface) = radius(on_interface).*sin(theta);
        end

        x(active) = x_old + dx;
        y(active) = y_old + dy;
        active = active(max(abs([x(active), y(active)]), [], 2) < 1-tol);
    end

    assert(all(isfinite(estimate)), '路径估计出现非有限数值。');
    total = cumsum(estimate);
    solution(k) = total(end)/num_path;
    solution_M(k,:) = (total(M_list)./M_list(:)).';
end
end

