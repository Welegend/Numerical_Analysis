function res=Romberg_calculus(fx,a,b)
% \param fx: syms expression
% \param a:  left  bound
% \param b:  right bound
% \return res: the root value of fx
%
% can change precision and max_iteration_times and min_diff_value at the following
precision=1e-6;
max_iteration_times=20;
min_diff_value=1e-4;
% ------------------------------------------------------------%
% example:
% fx=@(x)x^2*exp(x);
% a=0;
% b=0.8;
% Romberg_new(fx,a,b);
%
% return
% the answer is
%     0.5697
%     2.0890
%     0.2099
%     0.3196
%     0.3153
%     0.3148
%     0.3146
%     0.3146
% ------------------------------------------------------------%
eps=1e-6;
buf=zeros(10,10);
i=1;
h=(b-a);
buf(1,1)=(b-a)./2.*(fx(b)-fx(a));
buf(1,2)=h/2*(fx(b)+fx(a)+2*fx(a+h*i));
buf(2,1)=(4*buf(1,2)-buf(1,1))/(4-1);
if (abs((buf(i+1,1)-buf(i,1)))<eps)
    disp(buf(i+1,1));
else
    while (abs((buf(i+1,1)-buf(i,1)))>=eps)
        i=i+1;
        for m=2:i+1
            k=i-m;
            n=2.^i;
            h=(b-a)/2.^i;
            buf(1,i+1)=h/2*(fx(a)+fx(b));
            for j=1:n-1
                buf(1,i+1)=buf(1,i+1)+h/2*(2*fx(a+j*h)); % get T_1_i+1 value
            end
            buf(m+1,k+1)=(4^m*buf(m,k+2)-buf(m,k+1))/(4^m-1);
        end
    end
end
disp(buf);
res=buf(a+1,b+1);
end