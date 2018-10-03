% 函数功能：实现Ax=b的列主元高斯消去法（Gauss Elimination），其中A为n阶方阵
% 输入：矩阵A，b
% 输出：得到的解矩阵x
% 使用范围：列主元高斯消去法方法稳定，计算量小，适用于求解良态方程组（n<=30）
% 顺序主子式均不等于0，对称正定矩阵，严格对角占优矩阵时可以使用

function x = GE_equ(A, b)
%%  列主元高斯消去法，变换矩阵
Ab = [A, b];
[~, n] = size(A);
for k = 1: n - 1
     % 列主元，找到第k列最大值的那一行，并交换到对角线上
    [~, w] = max(abs(Ab(k: n, k)));
    Ab([k, w + k - 1], :) = Ab([w + k - 1, k], :);
    
     % 作用是用第k行的值往下消
    Ab(k + 1: n, k) = Ab(k + 1: n, k) / Ab(k, k);
    Ab(k + 1: n, k + 1: end) = Ab(k + 1: n, k + 1: end) - Ab(k + 1: n, k) * Ab(k, k + 1: end);
end

A = Ab(:, 1: size(A, 2));
b = Ab(:, end - size(b, 2) + 1: end);

%% 以下是A为方阵且上三角矩阵时，求解方程组Ax=b
x = UTri_equ(A, b);

end