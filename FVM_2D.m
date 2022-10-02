clc
close all
clear 

ax = 0; bx = 1;
ay = 0; by = 1;
vx = 1; vy = 1;

nx=20;ny=20;
dx=(bx-ax)/nx;
dy=(by-ay)/ny;
dt=0.004;
nt=1/dt;

x = ax+dx/2 : dx : bx;
y = ay+dy/2 : dy : by;
t=0:dt:1;
[X,Y] = meshgrid(x,y);

%Инициализация матрицы решений
u=zeros(nx,ny,nt);
u0=zeros(nx,ny);
u_analit=zeros(nx,ny,nt);

g=@(x,y,t) sin(x+y^2)+x*t;
f=@(x,y,t) x+t+cos(x+y^2)*(1+2*y);
bc_1=@(x,t)sin(x)+x*t;
bc_2=@(y,t)sin(y^2);
bc_3=@(x,t)sin(x+1)+x*t;
bc_4=@(y,t)sin(y^2+1)+t;
ic = @(x,y)sin(x+y^2);

%Начальные условия
for i = 1:(nx)
    for j = 1:(ny)
        u(i,j,1) = ic(x(i),y(j));
    end
end

%Граничные условия
for k=1:nt
    for i=1:nx
        u(i,1,k)=bc_1(x(i),t(k));   
    end
    for j=1:ny
        u(1,j,k)=bc_2(y(j),t(k));
    end
end

%Расчет и вывод аналитического решения
for k=1:nt
    for i = 1:nx
        for j = 1:ny
            u_analit(i,j,k) = g(x(i),y(j),t(k));
        end
    end
end
%Построение графика
subplot(2,1,1);
u_new_analit=u_analit(:,:,nt);
figure(1)
surf(X,Y,u_new_analit);
title('Аналитическоре решение')
%Создание меню
r = menu ( 'Menu',' Перенос  '  , ' Линейная интерполяция'  , ' Квадратичная интерполяция ');
 
switch r 
    
 case 1 
     %Перенос 
for k=1:nt
    for i=2:nx
        for j=2:ny  
            u(i,j,k+1)=u(i,j,k)-dt*vx*((u(i,j,k)-u(i-1,j,k))/dx)-vy*dt*((u(i,j,k)-u(i,j-1,k))/dy)+dt*f(x(i),y(j),t(k));  
        end
    end
end

    case 2
      
 for k=1:nt
    for i=1:nx
        u(i,ny,k)=bc_3(x(i),t(k));
    end
    for j=1:ny
        u(nx,j,k)=bc_4(y(j),t(k));
    end
end       
%Линейная интерполяция  
 for k=1:nt
    for i=2:nx-1
        for j=2:ny-1
            u(i,j,k+1)=u(i,j,k)-dt*vx*(u(i+1,j,k)-u(i-1,j,k))/(2*dx)-vy*dt*(u(i,j+1,k)-u(i,j-1,k))/(2*dy)+dt*f(x(i),y(j),t(k)); 
             
        end
    end
 end
       
    case 3
        
for k = 2:nt
    for i = 2:nx
    u(i,2,k) = g(x(i),y(2),t(k));
    u(i,3,k) = g(x(i),y(3),t(k));
    u(i,ny,k) = g(x(i),y(ny),t(k));
    end
    
    for j = 2:ny
    u(2,j,k) = g(x(2),y(j),t(k));
    u(3,j,k) = g(x(3),y(j),t(k));
    u(nx,j,k) = g(x(nx),y(j),t(k));
    end
end
%Квадратичная интерполяция
for k=1:nt-1
    for i=3:nx-1
        for j=3:ny-1
             u(i,j,k+1) = u(i,j,k)-dt*(((u(i-2,j,k)-7*u(i-1,j,k)+3*u(i,j,k)+3*u(i+1,j,k))*vx)/(dx*8) + ((u(i,j-2,k)-7*u(i,j-1,k)+3*u(i,j,k)+3*u(i,j+1,k))*vy)/(dy*8))+dt*f(x(i),y(j),t(k));
        end
    end
end

end

%Вывод решения
subplot(2,1,2);
u_new=u(:,:,nt);
% figure(2)
surf(X,Y,u_new);
title('Решение FVM')

% error=abs(u_analit(:,:,nt)-u(:,:,nt));
% figure(3)
% surf(X,Y,error);
% title('error')

