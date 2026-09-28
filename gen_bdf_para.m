
% GEN_BDF 生成 BDF-k 卷积求积系数和修正系数
% 保存为 .mat 文件供求解器使用
for idx_out = 1
% --- 参数设置 ---
    alpha_val = 0.05*idx_out;      % 分数阶导数阶数
    k_order = 4;          % BDF 阶数 (1-6)
    N_terms = 5000;       % 泰勒展开项数 (时间步数)
    
    fprintf('正在计算 BDF-%d (alpha=%.2f) 系数...\n', k_order, alpha_val);
    
    syms zeta z real
    
    
    delta_zeta = sym(0);
    for j = 1:k_order
        delta_zeta = delta_zeta + (1/j) * (1 - zeta)^j;
    end
    
    
    gen_func = delta_zeta^alpha_val;
    
    
    T = taylor(gen_func, zeta, 'Order', N_terms);
    
    
    cq_weights = zeros(N_terms, 1);
    
    [c_sym, t_sym] = coeffs(T, zeta);
    h_bar = waitbar(0);
    
    for i = 1:length(c_sym)
       
        p = double(polynomialDegree(t_sym(i), zeta));
        
        idx = p + 1;
        
        if idx <= N_terms
            cq_weights(idx) = double(c_sym(i));
        end
        if mod(i, 100) == 0
            waitbar(i / N_terms, h_bar, sprintf('已完成: %.1f%%', (i / N_terms) * 100));
        end
    end
    
    fprintf('   标准权重计算完成 (前5项): %s\n', mat2str(cq_weights(1:5), 4));
    
    % --- 2. 计算修正系数 (Correction Coefficients) ---
    % 求解 b_{j,m}
    correction_matrix = zeros(k_order-1, k_order-1);
    
    if k_order >= 2
        gamma_expr = zeta / (1 - zeta);
        
        for j = 0 : (k_order - 2)
            if j == 0
                curr_gamma = gamma_expr;
            else
                curr_gamma = zeta * diff(curr_gamma, zeta); 
            end
            
            % 目标函数: z^{-(j+1)} - gamma_j(e^{-z})/j!
            gamma_val = subs(curr_gamma, zeta, exp(-z));
            term1 = 1 / (z^(j+1));
            term2 = gamma_val / factorial(j);
            target_func = term1 - term2;
            
            target_series = taylor(target_func, z, 'Order', k_order);
            
            % 待定系数法
            syms c_var [1, k_order-1]
            ansatz = sym(0);
            for m = 1 : (k_order - 1)
                ansatz = ansatz + c_var(m) * exp(-m * z);
            end
            ansatz_series = taylor(ansatz, z, 'Order', k_order);
            
            eqns = [];
            for p = 0 : (k_order - 2)
                L = subs(diff(ansatz_series, z, p), z, 0);
                R = subs(diff(target_series, z, p), z, 0);
                eqns = [eqns, L == R];
            end
            
            sol = solve(eqns, c_var);
            sol_vec = struct2cell(sol);
            for m = 1 : (k_order - 1)
                correction_matrix(j+1, m) = double(sol_vec{m});
            end
        end
    end
    
    % --- 3. 保存 ---
    filename = sprintf('BDF%d_alpha%.2f_N%d.mat', k_order, alpha_val, N_terms);
    save(filename, 'alpha_val', 'k_order', 'cq_weights', 'correction_matrix');
    
    fprintf('成功保存系数到: %s\n', filename);
end