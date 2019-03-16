function res=newton_iteration(first_x_value,fx)
% \param x: start position
% \param fx: syms expression
% \return res: the root value of fx
%
% can change precision and max_iteration_times and min_diff_value at the following
precision=1e-6;
max_iteration_times=20;
min_diff_value=1e-4;
% ------------------------------------------------------------%
% example:
% >> syms x;
% fx=x-exp(-x);
% newton_iteration(0.5,fx);
%
% return
% the current x_value is 0.567143,
% the current iterate_times is 2
% ------------------------------------------------------------%

syms x;
diff_fx=diff(fx);
fx_next=x-fx/diff_fx;

x_before=double(first_x_value);
x_next=double(subs(fx_next,x,x_before));

iterate_times=1;
while (abs(x_next-x_before)>precision) && (iterate_times<=max_iteration_times)
    % if f'x < min_diff_value,
    % throw warning
    if abs(subs(diff_fx,x,x_before))<min_diff_value
        fprintf('\n diff less than min_diff_value \n min diff value is %06f\n', ...
            double(abs(subs(diff_fx,x,x_before))));
        break;
    end
    
    % print each iterate x_next value
    fprintf('\nthe current x_value is %06f,\nthe current iterate_times is %d\n', ...
        x_next,iterate_times);
    
    % next iterate
    x_before=x_next;
    x_next=double(subs(fx_next,x,x_before));
    iterate_times=iterate_times+1;
end
res=x_next;
end