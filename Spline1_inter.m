% 函数功能：三次样条插值函数，第一种边界条件：已知f''(a),f''(b)
% 输入：M0(即f''(a)),Mn（即f''(b)），n+1个插值点列向量x（从0到n）和y
% 输出：分段插值函数S(x)
function S = Spline1_inter(x, y, M0, Mn)
%% 三次样条插值的M值的方程组的导出
n = length(x);
u = (diff(x(1: n - 1))) ./ (diff(x(1: n - 1)) + diff(x(2: n))); % n-2维
la = 1 - u; % n-2维
f2 = (diff(y(2: n)) ./ diff(x(2: n)) - diff(y(1: n - 1)) ./ diff(x(1: n - 1))) ./ (diff(x(1: n - 1)) + diff(x(2: n))); % 二阶差商
d = 6 * f2;

%% 转化为三对角矩阵方程组
d(1) = d(1) - u(1) * M0;
d(end) = d(end) - la(end) * Mn;
M = Thomas_equ(u(2: end), 2 * ones(n - 2, 1), la(1: n - 3), d);

%% 代入S的方程组，S(x)为分段函数，段数不确定，用循环表示
M = [M0; M; Mn]; % M变成n维
S = cell(n - 1, 2); % 第一列存放插值表达式，第二列存放区间端点
for i = 1: n - 1
    %% 代入S的表达式
    S{i, 1} = @(xx)( ...
        (x(i + 1) - xx) .^3 ./ (6 .* (x(i + 1) - x(i))) .* M(i) ...
        + (xx - x(i)) .^3 ./ (6 .* (x(i + 1) - x(i))) .* M(i + 1) ...
        + (y(i) - (x(i + 1) - x(i)) .^2 / 6 .* M(i)) .* (x(i + 1) - xx) ./ (x(i + 1) - x(i)) ...
        + (y(i + 1) - (x(i + 1) - x(i)) .^2 ./ 6 .* M(i + 1)) .* (xx - x(i)) ./ (x(i + 1) - x(i)) ...
        ) .* (xx >= x(i) & xx <= x(i + 1)); % xx的定义域
    S{i, 2} = [x(i), x(i + 1)];
    
    %% 绘图
    fplot(S{i, 1}); % 函数句柄中的运算符、逻辑符需要转换成向量适用的(./ .* &)
    hold on
end

end