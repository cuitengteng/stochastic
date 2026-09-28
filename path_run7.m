function err = path_run7(alpha,beta,T,dt,num_path)%%%%
% 时空耦合：局部时间阶数为 1，非局部时间阶数为 beta，初值和外边界值均为零。
lx=-1; hx=1; ly=-1; hy=1;
num_mesh_x=20;
num_mesh_y=num_mesh_x;
refN=100;
x=linspace(lx,hx,num_mesh_x);
y=linspace(ly,hy,num_mesh_y);
[X,Y]=meshgrid(x,y);
x0=[X(:)';Y(:)'];
N=num_mesh_x*num_mesh_y;

alpha1=2;
% num_path=1600;
% alpha=0.;
% beta=0.5;
% dt=2e-4;
% T=1;
epoch=floor(num_path/dt);
assert(beta>0 && beta<=1,'要求 0<beta<=1。');
assert(alpha>0 && alpha<=2,'要求 0<alpha<=2。');
r=(dt*2^alpha*gamma(1+alpha/2)^2)^(1/alpha);
r1=(dt*2^alpha1*gamma(1+alpha1/2)^2)^(1/alpha1);

solution=zeros(N,1);
tic;
if alpha==2
    % WoS 传输：左侧 a=1，右侧 a=2，界面上解和通量连续。
    % dt 为平均出球时间上限；用平均值推进时间时钟属于离散近似。
    tol=1e-8;
    parfor k=1:N
        sol_stoc_new=0;
        for p=1:num_path
            x=x0(:,k);
            stime=0;
            estimate=0;
            for i=2:epoch
                if d_bound(x)<=tol || stime>=T
                    break;
                end
                if abs(x(1))<=tol
                    x(1)=0;
                    a=1.5;
                    on_interface=true;
                elseif x(1)<0
                    a=1;
                    on_interface=false;
                else
                    a=2;
                    on_interface=false;
                end
                radius=min(sqrt(4*a*dt),d_bound(x));
                if ~on_interface
                    radius=min(radius,abs(x(1)));
                end
                tau_q=radius^2/(4*a);
                if x(1)<=0
                    % 局部时钟按整数阶推进；界面点归入局部侧。
                    tau_q=min(tau_q,T-stime);
                    radius=sqrt(4*a*tau_q);
                    s_next=stime+tau_q;
                else
                    s_next=stime+rand_stable(beta,1,tau_q);
                end
                if s_next>T
                    break;
                end

                % 与 path_run6 相同的球内 Green 函数源项积分。
                r_source=radius*sqrt(rand*rand);
                theta_source=2*pi*rand;
                x_source=x+r_source*[cos(theta_source);sin(theta_source)];
                estimate=estimate+f_ins(x_source,alpha)*tau_q;

                if on_interface
                    theta_new=pi*(rand-0.5);
                    if rand<2/3
                        side=1;
                    else
                        side=-1;
                    end
                    x=x+radius*[side*cos(theta_new);sin(theta_new)];
                else
                    r_new=H(rand,alpha1,radius);
                    theta_new=2*pi*rand;
                    x=x+r_new*[cos(theta_new);sin(theta_new)];
                end
                stime=s_next;
                if i==epoch
                    error('WoS 路径达到步数上限，需增加 epoch。');
                end
            end
            sol_stoc_new=sol_stoc_new+estimate;
        end
        solution(k)=sol_stoc_new/num_path;
    end
else
parfor k=1:N
    sol_stoc_new=0;
    for p=1:num_path
        x=x0(:,k);
        stime=0;
        estimate=0;

        for i=2:epoch
            if d_bound(x)<=0
                break;
            end
            if d_boundl(x)>0 || (x(1)==0 && x(2)>-1 && x(2)<1)
                stime=stime+dt;
            elseif d_boundnl(x)>0
                stime=stime+rand_stable(beta,1,dt);
            else
                break;
            end
            if stime>T
                break;
            end

            if d_boundl(x)>0 || (x(1)==0 && x(2)>-1 && x(2)<1)
                % 局部区：布朗出球，跨界面反射，再尝试辅助跳跃。
                estimate=estimate+f_ins(x,alpha1)*dt;
                r_new=H(rand,alpha1,r1);
                theta_new=rand*2*pi;
                z=x+r_new.*[cos(theta_new);sin(theta_new)];

                if d_bound(z)<=0
                    break;
                elseif z(1)>=0 && z(1)<1 && z(2)>-1 && z(2)<1
                    z(1)=-abs(z(1));
                end

                r_new=H(rand,alpha,r);
                theta_new=rand*2*pi;
                y_cand=z+r_new.*[cos(theta_new);sin(theta_new)];
                if d_boundnl(y_cand)>0
                    x=y_cand;
                else
                    x=z;
                end

            elseif d_boundnl(x)>0
                % 非局部区：完整跳跃，再增加一次只接受区内落点的跳跃。
                estimate=estimate+f_ins(x,alpha)*dt;
                r_new=H(rand,alpha,r);
                theta_new=rand*2*pi;
                b=x+r_new.*[cos(theta_new);sin(theta_new)];

                if d_bound(b)<=0
                    break;
                elseif d_boundl(b)>0
                    x=b;
                elseif d_boundnl(b)>0
                    r_new=H(rand,alpha,r);
                    theta_new=rand*2*pi;
                    y_cand=b+r_new.*[cos(theta_new);sin(theta_new)];
                    if d_boundnl(y_cand)>0
                        x=y_cand;
                    else
                        x=b;
                    end
                else
                    break;
                end
            else
                break;
            end
        end
        sol_stoc_new=sol_stoc_new+estimate;
    end
    solution(k)=sol_stoc_new/num_path;
end

end

u_exact=main_st1(alpha,beta,T,1e-2,refN);
u_exact=u_exact';
solution_exact=u_exact(:);
err=calc_error_and_plot(solution,solution_exact,num_mesh_x,refN)
toc;


