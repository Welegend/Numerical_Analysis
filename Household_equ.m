% �������ܣ��ⷽ����Ax=b��Household QR�ֽ�ⷨ
% ���룺����A��b
% �����������x

function x = Household_equ(A, b)

[Q, R] = Household_QR([A, b]);
y = R(:, end);
R(:, end) = []; % ת��Ϊ���R * x = y��RΪ�����Ǿ���
end