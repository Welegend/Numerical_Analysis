% �������ܣ���ⷽ����Ax=b������ķ�����ΪAx=b
% ���룺����A��b��AΪn�������Ǿ���(Upper triangular matrix)��bΪ�о���
% �������x

function x = UTri_equ(A, b)
[~, n] = size(A);
x_temp = zeros(n, size(b, 2)); % Ԥ��x������
x = zeros(n, size(b, 2));

for k = 1: n
    x_temp(end, :) = b(end, :) ./ A(end, end); % ���㵱ǰ�����һ��x
    b = b - A(:, end) * x_temp(end, :); % �Ѽ��������x����뷽�̣��򻯷���
    x(n - k + 1, :) = x_temp(n - k + 1, :);
    
    % ÿ���һ��x�ͰѾ���ά
    x_temp(end, :) = [];
    A(:, end) = [];
    A(end, :) = [];
    b(end, :) = [];
end

end