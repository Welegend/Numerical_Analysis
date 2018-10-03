% 函数功能：求解方程组Ax=b，待解的方程组为Ax=b
% 输入：矩阵A，b；A为n阶下三角矩阵(Lower triangular matrix)，b为列矩阵
% 输出：解x

function x = LTri_equ(A, b)
[~, n] = size(A);
x_temp = zeros(n, size(b, 2)); % 预设x解向量
x = zeros(n, size(b, 2));

for k = 1: n
    x_temp(1, :) = b(1, :) ./ A(1, 1); % 先算当前的第1个x
    b = b - A(:, 1) * x_temp(1, :); % 把计算出来的x解代入方程，简化方程
    x(k, :) = x_temp(1, :);
    
    % 每解出一个x就把矩阵降维
    x_temp(1, :) = [];
    A(:, 1) = [];
    A(1, :) = [];
    b(1, :) = [];
end

end