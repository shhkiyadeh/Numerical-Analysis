% Sanaz Hami Hassan Kiyadeh, Test 2, Question 4:
clc
clear
a=-50;
b=50;
dt=0.01;
tmax=17;
c=1;
h=0.05;
x=a-h:h:b+h;
N=length(x)-2;
t=0;
u1=exp(-20*(x-2).^2)+exp(-(x-5).^2);
u2=u1;
w1=u1;
w2=u1;
M=tmax/dt;
% Upwind
for i=1:M
    for j=2:N+2
        w1(j)=u1(j)-c*(dt/h)*(u1(j)-u1(j-1));
    end
    for j=2:N+1
        w2(j)=0.5*(u2(j+1)+u2(j-1))+(c*dt/(2*h))*(u2(j+1)-u2(j-1));
    end
    t=t+dt;
    u1=w1;
    u2=w2;
end
% Lax-Fredrich
exact=exp(-20*(x-2-c*tmax).^2)+exp(-(x-5-c*tmax).^2);
figure(1)
plot(x,exact);
hold on
plot(x,u1,'o');
plot(x,u2,'*');
xlabel('x')
ylabel('u(x,17)')
legend('exact','Upwind method','Lax-Fredrich method')
title('Upwind & Lax-Fredrich method')

figure(2)
error1=(exact-u1).^2;
plot(x,error1)
xlabel('x')
ylabel('error')
title('Upwind method')

figure(3)
error2=(exact-u2).^2;
plot(x,error2)
xlabel('x')
ylabel('error')
title('Lax-Fredrich method')