% �������ܣ��ⷽ����Ax=b��Household QR�ֽ�ⷨ
% ���룺����A��b
% �����������x

function x = Household_equ(A, b)

[Q, R] = Household_QR([A, b]);
y = R(:, end);
R(:, end) = []; % ת��Ϊ���R * x = y��RΪ�����Ǿ���

[m, n] = size(R);
R = R(1: n, :); % ��m>n������£���R�����ɷ���
y = y(1: n); % R��ɷ���y��֮�仯
x_temp = zeros(n, 1); % Ԥ��x������
x = zeros(n, 1);

for k = 1: n
    x_temp(end) = y(end) / R(end, end); % �����n��x
    y = y - R(:, end) * x_temp(end); % �Ѽ��������x����뷽�̣��򻯷���
    x(n - k + 1) = x_temp(n - k + 1);
    x_temp(end) = [];
    R(:, end) = [];
    R(end, :) = [];
    y(end) = [];
end

end