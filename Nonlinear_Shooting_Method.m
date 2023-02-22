% Sanaz Hami Hassan Kiyadeh, Test 2, Question 1:

clc
clear
close all
a=0;
b=pi/2;
alpha=1;
beta=exp(1);
TOL=1e-4;
N=30;
h=(b-a)/N;
K=1;
M=100;
TK=(beta-alpha)/(b-a);
f=@(x,y,yp)yp.*cos(x)-y*log(y);
fy=@(y)-log(y)-1;
fyp=@(x)cos(x);

while K <= M
    w(1,1)=alpha;
    w(2,1)=TK;
    u1=0;
    u2=1;
    for i=1:N
        x=a+(i-1)*h;
        k(1,1)=h*w(2,i);
        k(1,2)=h*f(x,w(1,i),w(2,i));
        k(2,1)=h*(w(2,i)+1/2*k(1,2));
        k(2,2)=h*f(x+h/2,w(1,i)+1/2*k(1,1),w(2,i)+1/2*k(1,2));
        k(3,1)=h*(w(2,i)+1/2*k(2,2));
        k(3,2)=h*f(x+h/2,w(1,i)+1/2*k(2,1),w(2,i)+1/2*k(2,2));
        k(4,1)=h*(w(2,i)+k(3,2));
        k(4,2)=h*f(x+h,w(1,i)+k(3,1),w(2,i)+k(3,2));
        w(1,i+1)=w(1,i)+(k(1,1)+2*k(2,1)+2*k(3,1)+k(4,1))/6;
        w(2,i+1)=w(2,i)+(k(1,2)+2*k(2,2)+2*k(3,2)+k(4,2))/6;
        kp(1,1)=h*u2;
        kp(1,2)=h*fy(w(1,i))*u1+fyp(x)*u2;
        kp(2,1)=h*(u2+(1/2)*kp(1,2));
        kp(2,2)=h*(fy(w(1,i))*(u1+(1/2)*kp(1,1))+fyp(x+h/2)*(u2+(1/2)*kp(1,2)));
        kp(3,1)=h*(u2+(1/2)*kp(2,2));
        kp(3,2)=h*(fy(w(1,i))*(u1+(1/2)*kp(2,1))+fyp(x+h/2)*(u2+(1/2)*kp(2,2)));
        kp(4,1)=h*(u2+kp(3,2));
        kp(4,2)=h*(fy(w(1,i))*(u1+kp(3,1))+fyp(x+h)*(u2+kp(3,2)));
        u1=u1+(1/6)*(kp(1,1)+2*kp(1,2)+2*kp(3,1)+kp(4,1));
        u2=u2+(1/6)*(kp(1,2)+2*kp(2,2)+2*kp(3,2)+kp(4,2));
    end
    if abs(w(1,end)-beta)<=TOL
        for i=0:N
            x=a+i*h;
        end
    end
    TK=TK-((w(1,end)-beta)/u1);
    K=K+1; 
end
approximate=w(1,:);
xi=(0:h:pi/2);
exact=exp(sin(xi));
error=sum((approximate-exact).^2); 
e=(approximate-exact).^2;

figure(1)
plot(xi,approximate,'-');
hold on
plot(xi,exact,'*');
legend('approximate','exact');
xlabel('x');
ylabel('y(x)');
title('Nonlinear Shooting Method');

figure(2)
plot(xi,e)
xlabel('x');
ylabel('error')
title('Nonlinear Shooting Method');



