function res=Runge_Kutta(fx,a,b,alpha)
% \param fx: syms expression
% \param a:  left  bound
% \param b:  right bound
% \param alpha:  y0 value
% \return res: y_value list
%
% can change iterate times at the following
N=10;
x_list=zeros(1,15);   % size
y_list=zeros(1,15);   % size
% ------------------------------------------------------------%
% example:
% syms x y;
% fx=y;
% a=0;
% b=1;
% alpha=1;
% Runge_Kutta(fx,a,b,alpha);
%
% return
% the answer is
%          0         0.1000    0.2000    0.3000    0.4000    0.5000  
%          1.0000    1.1052    1.2214    1.3499    1.4918    1.6487

%          0.6000    0.7000    0.8000    0.9000    1.0000         0         
%          1.8221    2.0138    2.2255    2.4596    2.7183         0         
        
% ------------------------------------------------------------%

h=(b-a)/N;
syms x y;

k1=h*fx;
k2=simplify(h*subs(subs(fx,x,x+h/2),y,y+k1/2));
k3=simplify(h*subs(subs(fx,x,x+h/2),y,y+k2/2));
k4=simplify(h*subs(subs(fx,x,x+h),y,y+k3));

x_list(1)=a;
y_list(1)=alpha;

coefficient=simplify(y+1/6*(k1+2*k2+2*k3+k4));
disp([k1,k2,k3,k4,coefficient]);
for i=2:N+1
    x_list(i)=x_list(1)+(i-1)*h;
    y_list(i)=subs(subs(coefficient,x,x_list(i)),y,y_list(i-1));
end
disp([x_list;y_list]);
res= y_list;
end

    