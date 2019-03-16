function res=Legendre()
% P1(X)=1
% P2(X)=X
% ......
% pN+2(X)= (2*n + 3) / (n + 2) * x * PN+1(x) - (n + 1) / (n + 2) * PN(x)
%
% can change precision and max_iteration_times and min_diff_value at the following
max_iteration_times=20;
% n=1...max_iteration_times
% ------------------------------------------------------------%
% example:
% >> syms x;
% fx=x-exp(-x);
% newton_iteration(0.5,fx);
%
% return
% the current x_value is 0.567143,
% the current iterate_times is 2
% --

syms x;
coefficient_list=zeros(0,'syms',max_iteration_times);
coefficient_list(1)=1;
coefficient_list(2)=x;

for n=1:max_iteration_times-2
    coefficient_list(n+2)= simplify(...
        (2 * n + 3) / (n + 2) * x * coefficient_list(n+1) ...
        - (n + 1) / (n + 2) * coefficient_list(n)   ...
        );
end
pass=1
end