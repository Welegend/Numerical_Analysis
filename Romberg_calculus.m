function res=Romberg_calculus(fx,a,b)
% \param fx: syms expression
% \param a:  left  bound
% \param b:  right bound
% \return res: the fx calculus
% can change precision at the following
eps=1e-6;
buf=zeros(15);   % size
% ------------------------------------------------------------%
% example:
% fx=@(x)x^2*exp(x);
% a=0;
% b=0.8;
% Romberg_calculus(fx,a,b);
%
% return
% the answer is
%     0.5697
%     0.3172
%     0.3146
%     0.3146
%     0.3146
% ------------------------------------------------------------%
i=1;
buf(1,1)=(b-a)/2*(fx(b)+fx(a));
while (1)
    i=i+1;
    buf(1,i)=1/2*buf(1,i-1);
    for j=1:2^(i-2)
        buf(1,i)=buf(1,i)+(b-a)/(2^(i-1))*fx(a+(2*j-1)*(b-a)/2^(i-1)); % get T_1_i+1 value
    end
    
    for m=2:i
        k=i+1-m;
        buf(m,k)=(4^(m-1)*buf(m-1,k+1)-buf(m-1,k))/(4^(m-1)-1);
    end
    if (abs((buf(i,1)-buf(i-1,1)))<eps) 
        break;
    end
end
fprintf('%.6f\n',buf(:,1));
res=buf(i,1);
end