% 函数功能：解方程组Ax=b，Household QR分解解法，此方程解法只适用于A列满秩的情况
% 输入：矩阵A、b
% 输出：解向量x
% 使用范围：Household变换的QR分解计算方程组十分稳定，常用于计算病态方程组，但计算量约是高斯消去法的2倍

function x = Household_equ(A, b)

[~, R] = Household_QR([A, b]);
y = R(:, end - size(b, 2) + 1: end); % 考虑到b可能不止一列的情况
R(:, end - size(b, 2) + 1: end) = []; % 转化为求解R * x = y，R为上三角矩阵

[~, n] = size(R);
R = R(1: n, :); % 在m>n的情况下，把R矩阵变成方阵
y = y(1: n, :); % R变成方阵，y随之变化

%% 以下是R为方阵且上三角矩阵时，求解方程组
x_temp = zeros(n, size(y, 2)); % 预设x解向量
x = zeros(n, size(y, 2));

for k = 1: n
    x_temp(end, :) = y(end, :) ./ R(end, end); % 先算第n个x
    y = y - R(:, end) * x_temp(end, :); % 把计算出来的x解代入方程，简化方程
    x(n - k + 1, :) = x_temp(n - k + 1, :);
    x_temp(end, :) = [];
    R(:, end) = [];
    R(end, :) = [];
    y(end, :) = [];
end

end