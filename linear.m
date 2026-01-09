function [A,B] = linear(q0, dq0, tau0, M, C, G)
    n = 6;

    % 定义符号变量
    syms q [n 1] real
    syms dq [n 1] real
    syms tau [n 1] real

    % 状态向量 x = [q; dq]
    x = [q; dq];

    % 非线性动态方程表达式
    ddq = (M^-1) * (tau - C*dq - G);
    f = [dq; ddq];   % 状态导数 f(x, u)

    % 计算雅可比
    A_sym = jacobian(f, x);
    B_sym = jacobian(f, tau);

    % 替换为平衡点
    A = double(subs(A_sym, [q; dq; tau], [q0; dq0; tau0]));
    B = double(subs(B_sym, [q; dq; tau], [q0; dq0; tau0]));

end
