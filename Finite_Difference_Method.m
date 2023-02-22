% Sanaz Hami Hassan Kiyadeh, Test 2, Question 3:

clc
clear
close all
a=0;
b=pi;
c=0;
d=pi/2;
m=100;
n=100;
TOL=0.004;
N=1000;
h=pi/100;
k=pi/200;
f=@(x,y)-2*(cos(x+y)+cos(x-y));
g=@(x,y)cos(x)*cos(y);
for i=1:n-1
    x(i)=a+i*h;
end
for j=1:m-1
    y(j)=c+j*k;
end
for i=1:n-1
    for j=1:m-1
        w(i,j)=0;
    end
end
lambda=h^2/k^2;
miu=2*(1+lambda);
l=1;
while l <= N
    z=(-h^2*f(x(1),y(m-1))+g(a,y(m-1))+lambda*g(x(1),d)+lambda*w(1,m-2)+w(2,m-1))/miu;
    NORM=abs(z-w(1,m-1));
    w(1,m-1)=z;
    for i=2:n-2
        z=(-h^2*f(x(i),y(m-1))+lambda*g(x(i),d)+w(i-1,m-1)+w(i+1,m-1)+lambda*w(i,m-2))/miu;
        if abs(w(i,m-1)-z)>NORM; NORM=abs(w(i,m-1)-z);
        end
        w(i,m-1)=z;
    end
    z=(-h^2*f(x(n-1),y(m-1))+g(b,y(m-1))+lambda*g(x(n-1),d)+w(n-2,m-1)+lambda*w(n-1,m-2))/miu;
    if abs(w(n-1,m-1)-z)>NORM; NORM=abs(w(n-1,m-1)-z);
    end
    w(n-1,m-1)=z;
    for j=m-2:-1:2
        z=(-h^2*f(x(1),y(j))+g(a,y(j))+lambda*w(1,j+1)+lambda*w(1,j-1)+w(2,j))/miu;
        if abs(w(1,j)-z)>NORM; NORM=abs(w(1,j)-z);
        end
        w(1,j)=z;
        for i=2:n-2
            z=(-h^2*f(x(i),y(j))+w(i-1,j)+lambda*w(i,j+1)+w(i+1,j)+lambda*w(i,j-1))/miu;
            if abs(w(i,j)-z)>NORM; NORM=abs(w(i,j)-z);
            end
            w(i,j)=z;
        end
        z=(-h^2*f(x(n-1),y(j))+g(b,y(j))+w(n-2,j)+lambda*w(n-1,j+1)+lambda*w(n-1,j-1))/miu;
        if abs(w(n-1,j)-z)>NORM; NORM=abs(w(n-1,j)-z);
        end
        w(n-1,j)=z;
    end
    z=(-h^2*f(x(1),y(1))+g(a,y(1))+lambda*g(x(1),c)+lambda*(w(1,2))+w(2,1))/miu;
    if abs(w(1,1)-z)>NORM; NORM=abs(w(1,1)-z);
    end
    w(1,1)=z;
    for i=2:n-2
        z=(-h^2*f(x(i),y(1))+lambda*g(x(i),c)+w(i-1,1)+lambda*w(i,2)+w(i+1,1))/miu;
        if abs(w(i,1)-z)>NORM; NORM=abs(w(i,1)-z);
        end
        w(i,1)=z;
    end
    z=(-h^2*f(x(n-1),y(1))+g(b,y(1))+lambda*g(x(n-1),c)+w(n-2,1)+lambda*w(n-1,2))/miu;
    if abs(w(n-1,1)-z)>NORM; NORM=abs(w(n-1,1)-z);
    end
    w(n-1,1)=z;
    if NORM <= TOL
        for i=1:n-1
            for j=1:m-1
                xx(i)=x(i);
                yy(j)=y(j);
                ww(i,j)=w(i,j);
            end
        end
    end
    l=l+1;
end
for i=1:99
    for j=1:99
        exact(i,j)=cos(x(i))*cos(y(j));
    end
end
approximate=w;
error1=(approximate-exact).^2; 
error2=abs(approximate-exact);

figure(1)
mesh(x,y,exact)
title('exact')
xlabel('x')
ylabel('y')
zlabel('exact u(x,y)')
title('Finite Difference method')
zlim([-1,max(max(exact))])

figure(2)
mesh(x,y,approximate)
xlabel('x')
ylabel('y')
zlabel('approximate u(x,y)')
title('Finite Difference method')
zlim([-1,max(max(approximate))])

figure(3)
mesh(x,y,error1)
xlabel('x')
ylabel('y')
zlabel('L_2 error')
title('Finite Difference method')

figure(4)
mesh(x,y,error2)
xlabel('x')
ylabel('y')
zlabel('absolute error')
title('Finite Difference method')