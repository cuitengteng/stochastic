function [u_final,info]=main_st1(alpha,beta,T_end,tau,N)
%MAIN_ST1 时空局部--非局部耦合问题的 BDF4 参考解。
% 空间矩阵直接使用 solve_coupling_fixed 的内部点矩阵。
% 局部时间阶数为 1，非局部时间阶数为 beta，初值及外边界值均为零。
% 输出为 ndgrid 排列，与原 main_st1 的输出顺序一致。

assert(beta>0 && beta<=1 && T_end>=0 && tau>0);
if alpha==2
    % -div(a grad u)：左侧 a=1，右侧 a=2。
    % 界面满足 u_l=u_nl、du_l/dn=2*du_nl/dn；外边界为零。
    % 共用稳态传输矩阵，界面面系数为调和平均值 4/3。
    [~,A,space]=solve_coupling_fixed(2,N,4096);
else
    % 0<alpha<2：局部反射 Laplacian+Q，非局部 F+Q。
    [~,A,space]=solve_coupling_fixed(alpha,N,4096);
end
n=N-2;
total=n*n;
u_final=zeros(N,N);
info=struct('space',space,'tau',tau,'Nt',0);
if T_end==0
    return;
end

% 时间网格和 BDF4 权重。
Nt=ceil(T_end/tau);
tau=T_end/Nt;
bdf_coeffs=[25/12;-4;3;-4/3;1/4];
if beta==1
    cq_weights=zeros(Nt+1,1);
    count=min(5,Nt+1);
    cq_weights(1:count)=bdf_coeffs(1:count);
    correction_matrix=zeros(3,3);
    correction_matrix(1,:)=[31/24,-7/6,3/8];
else
    filename=sprintf('BDF4_alpha%.1f_N5000.mat',beta);
    coefficient_file=which(filename);
    assert(~isempty(coefficient_file),'缺少 BDF4 系数文件。');
    data=load(coefficient_file,'cq_weights','correction_matrix');
    cq_weights=data.cq_weights(:);
    correction_matrix=data.correction_matrix;
    assert(numel(cq_weights)>=Nt+1,'BDF4 权重长度不足。');
    assert(abs(cq_weights(1)-(25/12)^beta)<1e-10,'BDF 权重与 beta 不匹配。');
end

% 内部点右端项；与路径使用同一个 f_ins。
X=space.X;
Y=space.Y;
source=zeros(total,1);
for row=1:total
    source(row)=f_ins([X(row);Y(row)],alpha);
end

% 不同空间行使用不同的时间算子。
local_mask=X(:)<0;
factor=zeros(total,1);
for row=1:total
    if local_mask(row)
        factor(row)=bdf_coeffs(1)/tau;
    else
        factor(row)=cq_weights(1)/tau^beta;
    end
end
LHS=A+spdiags(factor,0,total,total);
[L_mat,U_mat,P_mat]=lu(LHS);
u_history=zeros(total,Nt+1);

% 全历史卷积，前 3 步加入启动修正。
for k=1:Nt
    hist_sum=zeros(total,1);
    for j=1:k
        if j<=4
            hist_sum(local_mask)=hist_sum(local_mask) ...
                +bdf_coeffs(j+1)*u_history(local_mask,k-j+1)/tau;
        end
        hist_sum(~local_mask)=hist_sum(~local_mask) ...
            +cq_weights(j+1)*u_history(~local_mask,k-j+1)/tau^beta;
    end
    rhs=source-hist_sum;
    if k<4
        rhs=rhs+correction_matrix(1,k)*source;
    end
    u_history(:,k+1)=U_mat\(L_mat\(P_mat*rhs));
end

u_final(2:end-1,2:end-1)=reshape(u_history(:,end),n,n);
info.tau=tau;
info.Nt=Nt;
info.beta_local=1;
info.beta_nonlocal=beta;
info.local_mask=local_mask;
end
