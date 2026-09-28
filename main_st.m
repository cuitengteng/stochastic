function [u_final,info]=main_st(alpha,beta,T_end,tau,N)
%MAIN_ST 时空局部--非局部耦合问题的 BDF4 参考解。
% 空间矩阵直接使用 solve_coupling_fixed 的内部点矩阵。
% 两侧时间阶数均为 beta，初值及外边界值均为零。
% 输出为 ndgrid 排列，与原 main_st 的输出顺序一致。

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
filename=sprintf('BDF4_alpha%g_N5000.mat',beta);
coefficient_file=which(filename);
assert(~isempty(coefficient_file),'缺少 BDF4 系数文件。');
data=load(coefficient_file,'cq_weights','correction_matrix');
cq_weights=data.cq_weights(:);
correction_matrix=data.correction_matrix;
assert(numel(cq_weights)>=Nt+1,'BDF4 权重长度不足。');

% 内部点右端项；与路径使用同一个 f_ins。
X=space.X;
Y=space.Y;
source=zeros(total,1);
for row=1:total
    source(row)=f_ins([X(row);Y(row)],alpha);
end

factor=cq_weights(1)/tau^beta;
LHS=A+factor*speye(total);
[L_mat,U_mat,P_mat]=lu(LHS);
u_history=zeros(total,Nt+1);

% 全历史卷积，前 3 步加入启动修正。
for k=1:Nt
    hist_sum=zeros(total,1);
    for j=1:k
        hist_sum=hist_sum+cq_weights(j+1)*u_history(:,k-j+1);
    end
    rhs=source-hist_sum/tau^beta;
    if k<4
        rhs=rhs+correction_matrix(1,k)*source;
    end
    u_history(:,k+1)=U_mat\(L_mat\(P_mat*rhs));
end

u_final(2:end-1,2:end-1)=reshape(u_history(:,end),n,n);
info.tau=tau;
info.Nt=Nt;
end

