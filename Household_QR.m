% �������ܣ�ʵ��A�����Household QR�ֽ⣬����ķ�����ΪAx=b
% ���룺����A
% ������ֽ��ľ���Q��R

function [Q, R] = Household_QR(A)

[m, n] = size(A);
H_fin = eye(m); % Ϊ��֤H_fin�ܺ�A��ˣ�ά����m�����ڴ��ÿ��H���ܳ˻�
for k = 1: min(m - 1, n)
    a = A(k: m, k); % A�ĵ�k��
    sigma = - sign(A(k, k)) * norm(a, 2); % normΪ�����ķ���
    w = a - sigma * eye(m - k + 1, 1);
    alpha = sigma * (sigma - a(1));
    H = eye(m); % ���ڽ�ά������
    H(k: end, k: end) = eye(m - k + 1) - w * w' / alpha; % H = I - 2 * u * u' �ı�ʽ
    H_fin = H * H_fin; % ��ά�����׺��H�˵�һ��
    A = H * A;
end

Q = H_fin';
R = A;

% ����m>n���������Q R���������õ���Rʼ�պ�A�Ĵ�Сһ��
if m > n
    Q = Q(:, 1: n);
    R = R(1: n, :);
end

end