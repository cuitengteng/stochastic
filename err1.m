%% This script runs Monte carlo for error with various alpha or beta

% x: initial position
% alpha_vec: index alpha
% num_path: the number of paths
% epoch: the uumber of cycles of a single path
% N:the number of initial point

alpha=0.3;
%beta=0.4;
num_path=10000;
beta_vec=[1];
T=1;
%num_path = [50 100 200 400 800 1600];
h=[1e-1 5e-2 2e-2 1e-2 5e-3 2e-3];
%dt=1e-4;
error = zeros(3,6);
for i=1
    beta=beta_vec(i);
    for j=1:6   
    error(i,j) = path_run7(alpha,beta,T,h(j),num_path);
    %error(i,j) = path_run6(alpha,T,dt,num_path(j));
    end
end

% error=[0.5425    0.3137    0.1630    0.1015    0.0614    0.0348;
%     0.5000    0.2891    0.1488    0.0918    0.0558    0.0315;
%     0.4652    0.2582    0.1322    0.0819    0.0493    0.0276;];%%%alpha=1.6

convergence_rates = zeros(3);
for i = 1:3
    log_dt = log(h);
    log_err = log(error(i,:));
    p = polyfit(log_dt, log_err, 1);
    convergence_rates(i) = p(1);
end
%
figure;
loglog(h,error(1,:),'Color','r','Marker','o','markersize',8,...
    'markerfacecolor','w','LineWidth',2);
hold on
loglog(h,error(2,:),'Color','b','Marker','^','markersize',8,...
    'markerfacecolor','w','LineWidth',2);
hold on
loglog(h,error(3,:),'Color','m','Marker','d' ,'markersize',8,...
    'markerfacecolor','w','LineWidth',2);
% hold on
% loglog(h,error(4,:),'Color','g','Marker','v' ,'markersize',8,...
%     'markerfacecolor','w','LineWidth',2);
 hold on
loglog(h,0.6*h.^(1/2), '--k', 'LineWidth', 1.5, 'DisplayName', 'O(M^{-1/2})');

%handle = legend('$\alpha=0.3$','$\alpha=0.9$','$\alpha=2$','${M}^{-\frac{1}{2}}$');
handle = legend('$\beta=0.1$','$\beta=0.5$','$\beta=1$','$\mathcal{O}({\Delta t}^{1/2})$');
set(handle,'Interpreter','latex','FontSize',14,'Location','Northeast')
%set(gca,'FontSize',14);
%yticks(1e-4:1e-2:1e-2);
%xticks(0:0.01:0.1);
%   % 返回 [x_min, x_max]
x_range = xlim; xlim([x_range(1)*2, x_range(2)]);
y_range = ylim;ylim([y_range(1)*0.5, y_range(2)*1.5]);
grid on
ylabel('Error','FontSize',16);
xlabel('$\Delta t$','Interpreter','latex','FontSize',16);
title('$\alpha=0.3,M=10000$','Interpreter','latex', 'FontSize', 16);

% 
