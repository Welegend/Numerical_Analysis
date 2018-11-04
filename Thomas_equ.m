% 函数功能：采用追赶法求解三对角方程组Ax=d
% 输入：n-1维列向量c，n-1维列向量a，n维列向量b，n维列向量d
% 输出：n维列向量x
function x = Thomas_equ(a, b, c, d)

n = length(d);
for i = 1: n - 1
    a(i) = a(i) / b(i);
    b(i + 1) = b(i + 1) - a(i) * c(i);
    d(i + 1) = d(i + 1) - a(i) * d(i);
end

%% 求解具有上带宽2的方程组
x(n) = d(n) / b(n);
for i = n - 1: -1: 1
    x(i) = (d(i) - c(i) * x(i + 1)) / b(i);
end
x = x';
end