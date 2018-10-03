% 函数功能：实现Ax=b的改进平方根法，其中A为n阶对称正定矩阵，将其分解为LDL'的形式
% 输入：矩阵A，b
% 输出：得到的解矩阵x
% 使用范围：改进平方根法比平方根法（GG'）少了开方运算，比高斯消去法和Doolittle分解法（LU）少了一半的乘除法
% 是求解系数矩阵为对称正定矩阵的线性方程组最有效的方法之一

function x = LDL_equ(A, b)

[~, n] = size(A);
%% 方阵的LDL分解(改进平方根分解)
A(2: n, 1) = A(2: n, 1) / A(1, 1); % 先求第一列的L矩阵，第一行不变不用求
for k = 2: n % 第k行循环
    for j = 2: n % 第j列逐次计算
        if j < k % 此时计算L矩阵
            A(k, j) = A(j, k) / A(j, j); % 计算L矩阵，此时k>j，这一步是区分Doolittle分解的步骤，在对称正定矩阵中简化了
        else % 此时计算U矩阵
            A(k, j) = A(k, j) - A(k, 1: k - 1) * A(1: k - 1, j); % 计算U矩阵，此时k<=j
        end
    end
end

%% 以下转化为求LDL'x=b方程组，拆分求解
% 解下三角矩阵方程组Ly=b
L = eye(n);
for j = 2: n
    L(j, 1: j - 1) = A(j, 1: j - 1);
end
y = LTri_equ(L, b);
% 解对角线矩阵方程组Dz=y
z = y ./ diag(A);
% 解上三角矩阵方程组L'x=z
x = UTri_equ(L', z);

end