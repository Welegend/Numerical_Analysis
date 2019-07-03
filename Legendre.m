function res=Legendre(x_value_list,iterate_times)
% \param first_x_value: start position
% \param iterate_times: max N num of PN
% 
% \return res: the root value of Legendre
% ------------------------------------------------------------%
% P1(X)=1
% P2(X)=X
% ......
% pN+2(X)= (2*n + 3) / (n + 2) * x * PN+1(x) - (n + 1) / (n + 2) * PN(x)
% ------------------------------------------------------------%
% example:
% >> Legendre([0.93 -0.93 0.66 -0.66 0.23 -0.23],7);
%
% return 
% the answer is
% 0.932470
% -0.932470
% 0.661209
% -0.661209
% 0.238619
% -0.238619
% ------------------------------------------------------------%

syms x;
coefficient_list=sym(zeros(0,iterate_times));
coefficient_list(1)=1;
coefficient_list(2)=x;

for index=1:iterate_times-2
    n=index-1;
    coefficient_list(index+2)= simplify(...
        (2 * n + 3) / (n + 2) * x * coefficient_list(index+1) ...
        - (n + 1) / (n + 2) * coefficient_list(index)   ...
        );
end
disp(transpose(coefficient_list));
res=newton_iteration(x_value_list,coefficient_list(iterate_times));
end