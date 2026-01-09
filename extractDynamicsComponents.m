function [M, C, G] = extractDynamicsComponents(tau, q, dq, ddq)
% 提取拉格朗日动力学中的 M(q), C(q,dq)*dq, G(q)
%
% 输入：
%   tau  : 拉格朗日广义力（符号向量） τ = M(q)*ddq + C(q,dq)*dq + G(q)
%   q    : 广义坐标符号向量 [q1; q2; ...; qn]
%   dq   : 一阶导 [dq1; dq2; ...; dqn]
%   ddq  : 二阶导 [ddq1; ddq2; ...; ddqn]
%
% 输出：
%   M    : 质量矩阵 M(q)
%   C_dq : 科氏/离心项 C(q,dq)*dq
%   G    : 重力项 G(q)

    n = length(q);
    M = sym(zeros(n));
    
    % 提取质量矩阵 M(q)
    for i = 1:n
        for j = 1:n
            M(i,j) = diff(tau(i), ddq(j));
        end
    end
    M = simplify(M);
    
    % 提取重力项 G(q)
    G = subs(tau, [dq; ddq], zeros(2*n, 1));
    G = simplify(G);
    
    % 提取科氏/离心项 C(q,dq)*dq
    C_dq = simplify(tau - M * ddq - G);

    C = sym(zeros(n));
    for i = 1:n
        for j = 1:n
            % 提取 C_dq(i) 中 dq(j) 的系数
            coeffs_ij = coeffs(C_dq(i), dq(j));
            if length(coeffs_ij) == 2
                C(i,j) = coeffs_ij(2);  % 有 dq(j) 项
            elseif isscalar(coeffs_ij) && has(C_dq(i), dq(j))
                C(i,j) = coeffs_ij(1);  % 只有一项，且是 dq(j)
            else
                C(i,j) = 0;
            end
        end
    end

    C = simplify(C);
end
