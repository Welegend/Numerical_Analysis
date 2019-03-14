function res=newton_iteration(first_x_value,fx)
% \param x: start position
% \param fx: syms expression
% \return res: the value the root in 
% 
%   can change precision and max_iteration_times at the following
disp(fx);
precision=1e-6;
max_iteration_times=100000;
% ------------------------------------------------------------%
syms x;
diff_fx=diff(fx);
f_x_next=(x-fx/diff_fx);

x_before=first_x_value;
x_next=subs(f_x_next,x,first_x_value);

iterate_times=1;
while (abs(x_next-x_before)>precision) && (iterate_times<=max_iteration_times)
    current_display=[x_next;iterate_times];
    formatSpec = '\nthe current x_value is %06f,\nthe current iterate_times is %d';
    fprintf(formatSpec,current_display);
    x_before=x_next;
    x_next=double(subs(f_x_next,x,x_before));
    iterate_times=iterate_times+1;
end
res=x_next;
end