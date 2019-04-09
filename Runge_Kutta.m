function res=Runge_Kutta(fx,a,b,alpha,N)
h=(b-a)/N;
syms x y;
k1=h*fx;
k2=sliplify(h*subs(subs(fx,x,x+h/2),y,y+k1/2));
k3=sliplify(h*subs(subs(fx,x,x+h/2),y,y+k2/2));
k4=sliplify(h*subs(subs(fx,x,x+h),y,y+k3));

x_list=zeros(1,15);
y_list=zeros(1,15);
x_list(1)=a;
y_list(1)=b;

coefficient=simplify(y_list(1)+1/6*(k1+2*k2+2*k3+k4));
for i=2:N
    x_list(i)=x_list(1)+(i-1)*h;
    y_list(i)=subs(subs(coefficient,x,x_list(i)),y,y_list(i-1));
end
disp([x_list;y_list]);
res= y_list;
end

    