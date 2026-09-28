%% Hybrid trajectory: WoS in Omega_nl + continuous Brownian path in Omega_l
clear; close all; clc;

figure('Position', [100, 100, 800, 550]);
hold on; axis equal; axis off;

%% 1. Basic setting
fill([-2, 0, 0, -2], [-1.5, -1.5, 1.5, 1.5], [0.96, 0.98, 1], ...
    'EdgeColor', 'k', 'LineWidth', 0.8);
fill([0, 2, 2, 0], [-1.5, -1.5, 1.5, 1.5], [0.98, 1, 0.97], ...
    'EdgeColor', 'k', 'LineWidth', 0.8);

plot([0, 0], [-1.5, 1.5], 'r-', 'LineWidth', 2.2);

text(-1, 1.2, '$\Omega_{n\ell}$', ...
    'Interpreter', 'latex', 'FontSize', 20, 'Color', 'b');

text(1, 1.2, '$\Omega_{\ell}$', ...
    'Interpreter', 'latex', 'FontSize', 20, 'Color', [0, 0.5, 0]);

theta = linspace(0, 2*pi, 100);

%% 2. Omega_nl: WoS balls for alpha-stable jump process
% For alpha < 2, the exit location may lie outside the ball,
% so the next center does not have to be on the previous circle.
nl_centers = [-1.5, 0.6;
              -0.9, 0.1;
              -0.3, 0.5];

r_nl = 0.22;   % fixed WoS radius

for i = 1:3
    plot(nl_centers(i,1) + r_nl*cos(theta), ...
         nl_centers(i,2) + r_nl*sin(theta), ...
         'b-', 'LineWidth', 1.3);

    text(nl_centers(i,1), nl_centers(i,2), ...
         sprintf('$X^{\\alpha}(t_%d)$', i-1), ...
         'Interpreter', 'latex', 'FontSize', 12, 'Color', 'b', ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle');
end

% Jump skeleton in Omega_nl
for i = 1:2
    dx = nl_centers(i+1,1) - nl_centers(i,1);
    dy = nl_centers(i+1,2) - nl_centers(i,2);

    quiver(nl_centers(i,1), nl_centers(i,2), dx, dy, ...
           0, 'b', 'LineWidth', 1.8, 'MaxHeadSize', 0.7);
end

%% 3. Coupling transition across Gamma
% The particle jumps from Omega_nl to Omega_l
x_cross_start = nl_centers(3,:);
x_local_start = [0.28, 0.55];

dx_couple = x_local_start(1) - x_cross_start(1);
dy_couple = x_local_start(2) - x_cross_start(2);

quiver(x_cross_start(1), x_cross_start(2), dx_couple, dy_couple, ...
       0, 'r', 'LineWidth', 3, 'MaxHeadSize', 1.0);

%% 3. Ω_l: Brownian WoS skeleton for alpha = 2
% In the local region, the true Brownian path is continuous.
% The balls below represent only the Brownian WoS skeleton:
% X_{i+1} lies on the boundary of B_r(X_i).

l_r = 0.19;   % fixed Brownian WoS radius

% first local center after crossing Gamma
l_centers = zeros(4,2);
l_centers(1,:) = [0.25, 0.52];

% choose directions manually for a compact schematic trajectory
phi = [-0.55, -0.25, -0.45];

for i = 1:3
    l_centers(i+1,:) = l_centers(i,:) + l_r*[cos(phi(i)), sin(phi(i))];
end

% plot local WoS balls
for i = 1:size(l_centers,1)
    plot(l_centers(i,1) + l_r*cos(theta), ...
         l_centers(i,2) + l_r*sin(theta), ...
         '-', 'LineWidth', 1.3, 'Color', [0, 0.6, 0]);

    text(l_centers(i,1)+0.03, l_centers(i,2)+0.04, ...
         sprintf('$X^{2}(t_%d)$', i), ...
         'Interpreter', 'latex', 'FontSize', 12, ...
         'Color', [0, 0.5, 0], ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle');
end

% connect the WoS centers
for i = 1:size(l_centers,1)-1
    dx = l_centers(i+1,1) - l_centers(i,1);
    dy = l_centers(i+1,2) - l_centers(i,2);

    quiver(l_centers(i,1), l_centers(i,2), dx, dy, ...
           0, 'Color', [0, 0.6, 0], ...
           'LineWidth', 1.8, 'MaxHeadSize', 0.7);
end

% optional: add a thin curve to indicate the hidden continuous Brownian path
tloc  = 1:size(l_centers,1);
ttloc = linspace(1, size(l_centers,1), 150);

bx = interp1(tloc, l_centers(:,1), ttloc, 'pchip');
by = interp1(tloc, l_centers(:,2), ttloc, 'pchip');

plot(bx, by, '--', 'Color', [0, 0.45, 0], 'LineWidth', 1.0);

%% 5. Interface label
text(0.06, -1.38, '$\Gamma$', ...
     'Interpreter', 'latex', ...
     'FontSize', 16, ...
     'HorizontalAlignment', 'left');

%% 6. Set range
xlim([-2.0, 2.0]);
ylim([-1.5, 1.5]);

hold off;