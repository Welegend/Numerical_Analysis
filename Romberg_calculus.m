function res=Romberg_calculus(fx,a,b)
% \param fx: syms expression
% \param a:  left  bound
% \param b:  right bound
% \return res: the fx calculas
%
% can change precision and max_iteration_times and min_diff_value at the following
eps=1e-6;
buf=zeros(15);   % size
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
        n=2.^i;
        h=(b-a)/(2^i);
        buf(1,i+1)=h/2*(fx(a)+fx(b));
        for j=1:n-1
            buf(1,i+1)=buf(1,i+1)+h/2*(2*fx(a+j*h)); % get T_1_i+1 value
        end
        
        for m=2:i+1
            k=i+2-m;
            buf(m,k)=(4^m*buf(m-1,k+1)-buf(m-1,k))/(4^m-1);
        end
    end
end
disp(buf(:,1));
res=buf(i,1);
end