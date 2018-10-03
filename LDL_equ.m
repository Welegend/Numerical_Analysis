% 函数功能：实现Ax=b的改进平方根法，其中A为m阶对称正定矩阵，将其分解为LDL'的形式
% 输入：矩阵A，b
% 输出：得到的解矩阵x
% 使用范围：改进平方根法比平方根法（GG'）少了开方运算，比高斯消去法和Doolittle分解法（LU）少了一半的乘除法
% 是求解系数矩阵为对称正定矩阵的线性方程组最有效的方法之一

function x = LDL_equ(A, b)
Ab = [A, b];
[m, n] = size(Ab);

%% 方阵的LDL分解(改进平方根分解)
Ab(2: m, 1) = Ab(2: m, 1) / Ab(1, 1); % 先求第一列的L矩阵，第一行不变不用求
for k = 2: m % 第k行循环
    for j = 2: n % 第j列逐次计算
        if j >= k % 此时计算U矩阵
            Ab(k, j) = Ab(k, j) - Ab(k, 1: k - 1) * Ab(1: k - 1, j); % 计算U矩阵，此时k<=j
        else % 此时计算L矩阵
            Ab(k, j) = Ab(j, k) / Ab(j, j); % 计算L矩阵，此时k>j，这一步是区分Doolittle分解的步骤，在对称正定矩阵中简化了
        end
    end
end

%% 以下转化为求DL'x=b方程组，拆分成L'x=z依次求解（第一个拆分在上一个步骤中已求得）
b = Ab(:, end - (n - m) + 1, end);
z = b ./ diag(Ab); % z是DL'
R = eye(m);
for j = 1: m - 1
    R(j, j + 1: m) = Ab(j + 1: m, j); % R矩阵是L'
end

%% 以下是R为m阶方阵且上三角矩阵时，求解方程组Rx=z
x_temp = zeros(m, size(z, 2)); % 预设x解向量
x = zeros(m, size(z, 2));
for k = 1: m
    x_temp(end, :) = z(end, :) ./ A(end, end); % 先算第n个x
    z = z - A(:, end) * x_temp(end, :); % 把计算出来的x解代入方程，简化方程
    x(m - k + 1, :) = x_temp(m - k + 1, :);
    x_temp(end, :) = [];
    A(:, end) = [];
    A(end, :) = [];
    z(end, :) = [];
end

end