% 函数功能：求解方程组Ax=b，待解的方程组为Ax=b
% 输入：矩阵A，b；A为n阶上三角矩阵(Upper triangular matrix)，b为列矩阵
% 输出：解x

function x = UTri_equ(A, b)
[~, n] = size(A);
x_temp = zeros(n, size(b, 2)); % 预设x解向量
x = zeros(n, size(b, 2));

for k = 1: n
    x_temp(end, :) = b(end, :) ./ A(end, end); % 先算当前的最后一个x
    b = b - A(:, end) * x_temp(end, :); % 把计算出来的x解代入方程，简化方程
    x(n - k + 1, :) = x_temp(n - k + 1, :);
    
    % 每解出一个x就把矩阵降维
    x_temp(end, :) = [];
    A(:, end) = [];
    A(end, :) = [];
    b(end, :) = [];
end

end