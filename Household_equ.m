% 函数功能：解方程组Ax=b，Household QR分解解法
% 输入：矩阵A、b
% 输出：解向量x

function x = Household_equ(A, b)

[Q, R] = Household_QR([A, b]);
y = R(:, end);
R(:, end) = []; % 转化为求解R * x = y，R为上三角矩阵
end