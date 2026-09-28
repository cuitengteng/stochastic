function rel2_err= calc_error_and_plot(u_py, u_ref, N_py, N_ref)
    % 1. 确保参考解 u_ref 是正确的 2D 形状
    % 防止读入的是 flatten 后的向量
    u_py=reshape(u_py, [N_py, N_py]);
    u_py=u_py';
    %u_ref=solution_exact;
    %N_py=25;
    %N_ref=100;
    if isvector(u_ref)
        u_ref = reshape(u_ref, [N_ref, N_ref]);
    end

    % 2. 定义坐标向量 (Grid Vectors)
    %这是 griddedInterpolant 最喜欢的格式
    x_py = linspace(-1, 1, N_py)';
    y_py = linspace(-1, 1, N_py)';

    x_ref = linspace(-1, 1, N_ref)';
    y_ref = linspace(-1, 1, N_ref)';

    % 3. 生成 Python 解的目标网格 (用于绘图和后续误差计算)
    [X_target, Y_target] = ndgrid(x_py, y_py);

    % 4. 构建插值器 (基于 Reference 网格)
    % 关键点：因为 u_ref 对应的是 ndgrid 格式 (行对应 x_ref, 列对应 y_ref)
    % 所以我们传入元胞数组 {x_ref, y_ref}
    F = griddedInterpolant({x_ref, y_ref}, u_ref', 'linear');

    % 5. 执行插值
    % 直接将目标网格矩阵传入 F
    u_ref_interp = F(X_target, Y_target);

    % 6. 计算误差
    error_abs = abs(u_py - u_ref_interp);
    max_err = max(error_abs(:));
    rel2_err = norm(error_abs(:)) / norm(u_ref_interp(:));
    %l2_err=sqrt(1/(N_py.^2)*sum((error_abs(:)).^2));
    % 7. 打印统计
    fprintf('=== 误差统计 ===\n');
    fprintf('网格映射: %dx%d -> %dx%d\n', N_ref, N_ref, N_py, N_py);
    fprintf('L_inf Error (Max): %.4e\n', max_err);
    fprintf('Rel. L_2 Error   : %.4e\n', rel2_err);
    %fprintf(' L_2 Error   : %.4e\n', l2_err);
    % 8. 绘图
    % figure('Position', [100, 500, 1200, 350]);
    % 
    % subplot(1,3,1);
    % % 这里的转置 ' 是为了让 imagesc 正确显示 x轴水平，y轴垂直
    % % 或者使用 pcolor(X_target, Y_target, u_py); shading nearest;
    % pcolor(X_target, Y_target, u_py); shading interp;
    % colorbar; title([' Solution (' num2str(N_py) 'x' num2str(N_py) ')']);
    % xlabel('x'); ylabel('y');
    % 
    % subplot(1,3,2);
    % pcolor(X_target, Y_target, u_ref_interp); shading interp;
    % colorbar; title(['Ref solution (' num2str(N_ref) '->' num2str(N_py) ')']);
    % xlabel('x'); ylabel('y');
    % 
    % subplot(1,3,3);
    % pcolor(X_target, Y_target, error_abs); shading interp;
    % colorbar; title(['Diff (Max: ' num2str(max_err, '%.1e') ')']);
    % colormap('jet');
    % xlabel('x'); ylabel('y');
    % 
    % sgtitle('Error Analysis');
%end