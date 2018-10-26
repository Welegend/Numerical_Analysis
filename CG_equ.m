% 函数功能：共轭梯度法是一种迭代解方程组方法，解方程组Ax=b，A是对称正定矩阵
% 输入：对称正定矩阵A, 列向量b
% 输出：方程组的解向量x

function x = CG_equ(A, b)

n = size(A, 2); % A为n阶对称正定矩阵
x = zeros(n, 1); % x为迭代初始值
r = b; % r为残差，因为x初始值为0，这里省去Ax
d = r; % d为搜索方向，第一次的值取负梯度，即函数下降最快方向

k = 0; % k为迭代次数，最大为n
er = 10e-8; % er为迭代精度

while norm(r) > er && k < n
    alpha = r' * r / (d' * A * d); % alpha为最佳步长
    x = x + alpha * d;
    rr = b - A * x; % 迭代后的残差
    beta = rr' * rr / (r' * r);
    d = rr + beta * d;

    r = rr; % 用完第k次的r后，舍弃掉
    k = k + 1; % 标志第k次迭代结束
end

end