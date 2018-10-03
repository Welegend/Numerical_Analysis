% �������ܣ���ⷽ����Ax=b������ķ�����ΪAx=b
% ���룺����A��b��AΪn�������Ǿ���(Lower triangular matrix)��bΪ�о���
% �������x

function x = LTri_equ(A, b)
[~, n] = size(A);
x_temp = zeros(n, size(b, 2)); % Ԥ��x������
x = zeros(n, size(b, 2));

for k = 1: n
    x_temp(1, :) = b(1, :) ./ A(1, 1); % ���㵱ǰ�ĵ�1��x
    b = b - A(:, 1) * x_temp(1, :); % �Ѽ��������x����뷽�̣��򻯷���
    x(k, :) = x_temp(1, :);
    
    % ÿ���һ��x�ͰѾ���ά
    x_temp(1, :) = [];
    A(:, 1) = [];
    A(1, :) = [];
    b(1, :) = [];
end

end