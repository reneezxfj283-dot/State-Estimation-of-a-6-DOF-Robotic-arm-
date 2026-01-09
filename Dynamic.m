function dyn = Dynamic()
syms q1 q2 q3 q4 q5 q6 real
syms dq1 dq2 dq3 dq4 dq5 dq6 real
syms ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 real
global T_all

dyn = struct(); 
q  = [q1;  q2;  q3;  q4;  q5;  q6];
dq = [dq1; dq2; dq3; dq4; dq5; dq6];
ddq = [ddq1; ddq2; ddq3; ddq4; ddq5; ddq6];
n = 6; % DOF

DH;
rob_Para;

for l=1:2                               
% Extracting robot arm parameters from the parameter structure
m = rob_para(l).masses; % 各连杆质量mass [kg]  
r = rob_para(l).r;
% L = rob_para.lengths; % 各连杆长度 [m]
% I_diag = rob_Para.inertias; % 各连杆惯性矩阵对角线元素

%% Kinetic energy
K = 0; % 总动能初始化
T = cell(1, n); 
I_pseudo = cell(1, n);  % 伪惯性矩阵初始化
T0 = eye(4);
for i = 1:n
    I_pseudo{i} = m(i) * r{i} * r{i}';
    if i == 1
        T{i} = T0 * T_all(:,:,i);   % 第一个关节
    else
        T{i} = T{i-1} * T_all(:,:,i);  % 后续关节
    end                           

    for j = 1:i
        dT_dqj = diff(T{i}, q(j));             %Ti对qj求偏导  
        for k = 1:i
            dT_dqk = diff(T{i}, q(k));         %Ti对qk求偏导
            dK_term = trace( dT_dqj(1:3,1:3) * I_pseudo{i} * dT_dqk(1:3,1:3).' );
            K = K + 0.5 * dK_term * dq(j) * dq(k);
        end
    end
end

K = simplify(K);
disp('动能表达式 K(q,dq)：');
pretty(K)

%% potential energy
g_vec = [0; 0; -9.81];   % 定义重力向量（基坐标系下）
P = sym(0);
T = cell(1, n);     
% 定义每个连杆的质量和质心位置
for i = 1:n
    % 第i个连杆的齐次变换矩阵（从 base 到该连杆）
    if i == 1
        T{i} = T0 * T_all(:,:,i);   % 第一个关节
    else
        T{i} = T{i-1} * T_all(:,:,i);  % 后续关节
    end                           

    pi = T{i}(1:3,1:3) * r{i} + T{i}(1:3,4);     % pi = T_i * r_i   % 提取第i个质心在世界坐标系下的位置向量。   Ti(1:3,1:3)为旋转部分，Ti(1:3,4)平移部分
    Pi = -m(i) * g_vec.' * pi;               % 当前连杆的势能项：-m_i * g^T * p_i

    % 累加总势能
    P = P + Pi;
end

P = simplify(P);
disp('总势能表达式 P(q)：');
pretty(P)

%% 动力学方程
L = K - P;    % Lagrangian

%构造动力学方程左侧
tau = sym(zeros(n,1));  % 初始化关节力矩

for i = 1:n
    dL_dqi  = diff(L, q(i));
    dL_ddqi = diff(L, dq(i));

    % 对时间求导：d/dt(∂L/∂dq_i)
    d_dt_dL_ddqi = 0;
    for j = 1:n
        d_dt_dL_ddqi = d_dt_dL_ddqi + diff(dL_ddqi, q(j)) * dq(j) + diff(dL_ddqi, dq(j)) * ddq(j);
    end

    % 拉格朗日方程
    tau(i) = simplify(d_dt_dL_ddqi - dL_dqi);
end

% 分离M、C、G，从 tau = M(q)*ddq + C(q,dq) + G(q) 提取结构
[M, C, G] = extractDynamicsComponents(tau, q, dq, ddq);   

%% 线性化
% q0   = [0;  0;  2*pi/3;  -pi/4;  -pi/2;  0];
q0   = [0;  0;  2.09;  -0.79;  -1.57;  0];
dq0  = zeros(6,1); 
G0   = double(subs(G, q, q0));  
tau0 = G0;   % 输入力矩 = 重力补偿

[A,B] = linear(q0, dq0, tau0, M, C, G);
% lin(l).A=A;
% lin(l).B=B;
%% 离散化
d = size(A,1);
M = [A, B;zeros(size(B,2), d + size(B,2))];
M_d = expm(M * 0.002);
Ad = M_d(1:d, 1:d);
Bd = M_d(1:d, d+1:end);

% 保存结果作为函数
dyn(l).Ad = Ad;
dyn(l).Bd = Bd;

end
end