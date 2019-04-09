function res=Romberg_new(fx,a,b)
eps=1e-5;
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