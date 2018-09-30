% 函数功能：解方程组Ax=b，Household QR分解解法
% 输入：矩阵A、b
% 输出：解向量x

function x = Household_equ(A, b)

[Q, R] = Household_QR([A, b]);
y = R(:, end);
R(:, end) = []; % 转化为求解R * x = y，R为上三角矩阵

[m, n] = size(R);
R = R(1: n, :); % 在m>n的情况下，把R矩阵变成方阵
y = y(1: n); % R变成方阵，y随之变化
x_temp = zeros(n, 1); % 预设x解向量
x = zeros(n, 1);

for k = 1: n
    x_temp(end) = y(end) / R(end, end); % 先算第n个x
    y = y - R(:, end) * x_temp(end); % 把计算出来的x解代入方程，简化方程
    x(n - k + 1) = x_temp(n - k + 1);
    x_temp(end) = [];
    R(:, end) = [];
    R(end, :) = [];
    y(end) = [];
end

end