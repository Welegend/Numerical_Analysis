% 函数功能：实现A矩阵的Household QR分解，待解的方程组为Ax=b
% 输入：矩阵A
% 输出：分解后的矩阵Q、R

function [Q, R] = Household_QR(A)

[m, n] = size(A);
H_fin = eye(m); % 为保证H_fin能和A相乘，维数是m，用于存放每次H的总乘积
for k = 1: min(m - 1, n)
    a = A(k: m, k); % A的第k列
    sigma = - sign(A(k, k)) * norm(a, 2); % norm为向量的范数
    w = a - sigma * eye(m - k + 1, 1);
    alpha = sigma * (sigma - a(1));
    H = eye(m); % 用于降维、升阶
    H(k: end, k: end) = eye(m - k + 1) - w * w' / alpha; % H = I - 2 * u * u' 的变式
    H_fin = H * H_fin; % 降维、升阶后把H乘到一起
    A = H * A;
end

Q = H_fin';
R = A;

% 考虑m>n的情况，简化Q R矩阵，这样得到的R始终和A的大小一样
if m > n
    Q = Q(:, 1: n);
    R = R(1: n, :);
end

end