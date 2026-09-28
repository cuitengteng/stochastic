function [u_exact,A,info] = solve_coupling_fixed(alpha,N,fftN)
%SOLVE_COUPLING_FIXED  main6 稳态局部--非局部耦合问题的 FDM 参考解。
%
% 单文件版本，无需其他辅助函数。
%
% 调用：
%   u_exact = solve_coupling_fixed(alpha,N)
%   [u_exact,A,info] = solve_coupling_fixed(alpha,N,fftN)
%
% 输入：
%   alpha : 分数阶，主耦合模型要求 0<alpha<2。
%   N     : 每个方向包含边界在内的节点数，必须为偶数。
%           偶数 N 使 x1=0 位于 -h/2 和 h/2 之间的网格面上。
%   fftN  : 计算分数中心差分系数的 FFT 尺寸，默认 4096。
%
% 输出：
%   u_exact : N^2 x 1，包含零边界，按 meshgrid 展平顺序输出，
%             可直接与 path_run6_remark26_commented.m 的 X(:),Y(:) 配合。
%   A       : (N-2)^2 x (N-2)^2 的内部离散矩阵，按 ndgrid 顺序。
%   info    : 网格、内部解、右端项、残差、对称性和内存信息。
%
% 对 0<alpha<2，离散算子为
%   Omega_l  : -Delta_h^{reflected} + Q_h,
%   Omega_nl : F_h + Q_h.
% F_h 是零外延的全空间分数中心差分算子；Q_h 是只限制在
% Omega_nl 内的第二份相互作用。因此非局部区内、跨区和域外系数
% 分别为 2、1、1。
%
% 对 alpha<2：外边界为齐次零 Dirichlet；局部侧在 x1=0 为反射
% Neumann 条件。非零外部数据需要另行加入右端项。
%
% alpha=2 solves -div(a grad u)=f, a=1 on the left and a=2 on the right.
% Interface: continuous u and a_left*du_left/dn=a_right*du_right/dn.
% No additional zero-Neumann condition is imposed for alpha=2.
% Outer boundary: g=0 (nonzero g requires boundary contributions to rhs).

if nargin<3 || isempty(fftN), fftN=4096; end
assert(isscalar(alpha) && alpha>0 && alpha<=2,'要求 0<alpha<=2。');
assert(isscalar(N) && N>=6 && N==floor(N) && mod(N,2)==0, ...
    'N 必须是至少为 6 的偶数，并且 N 包含边界节点。');
assert(isscalar(fftN) && fftN==floor(fftN) && fftN>=2*N, ...
    'fftN 必须是整数且至少为 2*N。');

n=N-2;                         % 每个方向的内部节点数
h=2/(N-1);
axis_full=linspace(-1,1,N);
axis_interior=axis_full(2:end-1);
[X,Y]=ndgrid(axis_interior,axis_interior);
[I,J]=ndgrid(1:n,1:n);
I=I(:); J=J(:);

local_mask=X(:)<0;
nonlocal_indices=find(~local_mask); % 偶数 N 保证不存在 x1=0 内部节点
rhs=1+0.5*cos(X(:).*Y(:));

total=n*n;
matrix_memory_gib=8*double(total)^2/1024^3;
fprintf('Assembling %d-by-%d system (dense equivalent %.3f GiB).\n', ...
    total,total,matrix_memory_gib);
if alpha==2
    A=spalloc(total,total,5*total);
else
    A=zeros(total,total);
end

if alpha<2
    coefficients=fcd_coefficients(alpha,n,fftN)/h^alpha;
end

for row=1:total
    i=I(row); j=J(row);

    if alpha==2
        % Conservative transmission flux: a_left=1, a_right=2.
        % Equal half-cell distances give interface conductance 4/3.
        % Missing exterior neighbors have zero Dirichlet values; their
        % conductance remains in the diagonal. No reflecting interface here.
        a=1+double(~local_mask(row));
        left=a; right=a;
        if i==n/2, right=4/3; end
        if i==n/2+1, left=4/3; end
        A(row,row)=(left+right+2*a)/h^2;
        if i>1, A(row,row-1)=-left/h^2; end
        if i<n, A(row,row+1)=-right/h^2; end
        if j>1, A(row,row-n)=-a/h^2; end
        if j<n, A(row,row+n)=-a/h^2; end
    elseif local_mask(row)
        % 局部区：反射五点 Laplacian。
        % 最后一条局部网格线 i=n/2 位于 x1=-h/2。
        % 用 u_ghost(h/2)=u(-h/2) 消去 ghost 后，对角元从 4/h^2
        % 变为 3/h^2，并且不存在跨界面的 Brownian 邻接边。
        interface_face=(i==n/2);
        A(row,row)=(4-double(interface_face))/h^2;
        if i>1, A(row,row-1)=-1/h^2; end
        if i<n && ~interface_face, A(row,row+1)=-1/h^2; end
        if j>1, A(row,row-n)=-1/h^2; end
        if j<n, A(row,row+n)=-1/h^2; end

    else
        % 非局部区：完整零外延分数算子 F_h。
        lookup=sub2ind([n,n],abs(I-i)+1,abs(J-j)+1);
        A(row,:)=coefficients(lookup)';
    end

    if alpha<2
        % 所有行都加入只指向 Omega_nl 的 Q_h。
        % local 行：这给出跨区耦合；nonlocal 行：这是第二份区内强度。
        lookup=sub2ind([n,n], ...
            abs(I(nonlocal_indices)-i)+1,abs(J(nonlocal_indices)-j)+1);
        weights=-coefficients(lookup)';
        weights(nonlocal_indices==row)=0;

        tolerance=1e-12*max(1,max(abs(weights)));
        weights(weights<0 & weights>-tolerance)=0;
        assert(all(weights>=0), ...
            '分数差分出现负跳跃权；请增大 fftN。');

        A(row,row)=A(row,row)+sum(weights);
        A(row,nonlocal_indices)=A(row,nonlocal_indices)-weights;
    end
end

interior_solution=A\rhs;
relative_residual=norm(A*interior_solution-rhs)/norm(rhs);

% info.grid_ndgrid 的第一维对应 x1；u_exact 转成原随机代码的
% meshgrid/X(:),Y(:) 展平顺序。
grid_ndgrid=zeros(N,N);
grid_ndgrid(2:end-1,2:end-1)=reshape(interior_solution,n,n);
u_exact=reshape(grid_ndgrid',[],1);

info=struct();
info.axis=axis_full;
info.X=X;
info.Y=Y;
info.grid_ndgrid=grid_ndgrid;
info.interior_solution=interior_solution;
info.rhs=rhs;
info.h=h;
info.fftN=fftN;
info.residual=relative_residual;
info.symmetry=norm(A-A','fro')/norm(A,'fro');
info.matrix_memory_gib=matrix_memory_gib;
info.matrix_ordering='interior ndgrid; first coordinate varies fastest';
info.solution_ordering='full-grid meshgrid vector, compatible with X(:),Y(:)';
if alpha==2
    info.interface_discretization='harmonic transmission flux, continuous trace';
    info.diffusion_coefficients=[1 2];
    info.interface_conductance=4/3;
else
    info.interface_discretization='face-centred reflected ghost';
end

fprintf('relative residual = %.3e, symmetry defect = %.3e\n', ...
    info.residual,info.symmetry);
end


function coefficients=fcd_coefficients(alpha,n,fftN)
% 分数中心差分符号：
% [4 sin(k1/2)^2 + 4 sin(k2/2)^2]^(alpha/2).
frequency=2*pi*(0:fftN-1)/fftN;
symbol=(4*sin(frequency(:)/2).^2+4*sin(frequency/2).^2).^(alpha/2);
full_coefficients=real(ifft2(symbol));
coefficients=full_coefficients(1:n,1:n);

off_diagonal=coefficients;
off_diagonal(1,1)=0;
if max(off_diagonal(:))>1e-10
    warning('检测到正的 FCD 非对角系数；建议增大 fftN。');
end
end
