% �������ܣ�����׷�Ϸ�������ԽǷ�����Ax=d
% ���룺n-1ά������c��n-1ά������a��nά������b��nά������d
% �����nά������x
function x = Thomas_equ(a, b, c, d)

n = length(d);
for i = 1: n - 1
    a(i) = a(i) / b(i);
    b(i + 1) = b(i + 1) - a(i) * c(i);
    d(i + 1) = d(i + 1) - a(i) * d(i);
end

%% �������ϴ���2�ķ�����
x(n) = d(n) / b(n);
for i = n - 1: -1: 1
    x(i) = (d(i) - c(i) * x(i + 1)) / b(i);
end
x = x';
end