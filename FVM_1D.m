clc
close all
clear 

ax = 0; bx = 1;vx = 2;

nx=20;
dx=(bx-ax)/nx;
dt=0.004;
nt=1/dt;

x = ax+dx/2 : dx : bx;
t=0:dt:1;

%Инициализация матрицы решений
u=zeros(nx,nt);
u_analit=zeros(nx,nt);
error=zeros(nx,nt);
g=@(x,t) 5*x^2+4*t;
f=@(x,t) 4+20*x;
bc=@(t)4*t;
bc_2=@(t)5+4*t;
ic = @(x)5*x^2;


%Начальные условия
for i = 1:nx
   u(i,1) = ic(x(i));   
end

%Граничные условия
for k=1:nt
     u(1,k)=bc(t(k));       
end

%Расчет и вывод аналитического решения
for k=1:nt
    for i = 1:nx
        u_analit(i,k) = g(x(i),t(k));
    end
end

%Построение графика
figure(1)
subplot(2,1,1);
contourf(u_analit,100,'linecolor','non');
colormap(jet(256));
colorbar
caxis([-10,10])
title('Аналитическоре решение')

%Создание меню
r = menu ( 'Menu',' Перенос '  , ' Линейная интерполяция'  , ' Квадратичная интерполяция ');
 
switch r 
    
 case 1 
for k=1:nt
    for i=2:nx 
        u(i,k+1)=u(i,k)-dt*vx*((u(i,k)-u(i-1,k))/dx)+dt*f(x(i),t(k));  
    end
end


case 2        
for k=1:nt
    u(nx,k)=bc_2(t(k));
end
     
%Линейная интерполяция  
 for k=1:nt
     for i=2:nx-1
        u(i,k+1)=u(i,k)-dt*vx*(u(i+1,k)-u(i-1,k))/(2*dx)+dt*f(x(i),t(k));      
     end  
 end
  
  case 3       
for k = 2:nt
    u(2,k) = g(x(2),t(k));
    u(3,k) = g(x(3),t(k));
    u(nx,k) = g(x(nx),t(k));
end
%Квадратичная интерполяция
for k=1:nt-1
    for i=3:nx-1
         u(i,k+1) = u(i,k)-dt*(((u(i-2,k)-7*u(i-1,k)+3*u(i,k)+3*u(i+1,k))*vx)/(dx*8))+dt*f(x(i),t(k));
     end

end

end

subplot(2,1,2);
%figure(2)
contourf(u,100,'linecolor','non');
colormap(jet(256));
colorbar
caxis([-10,10])
title('Решение FVM')

%Вычисление ошибки

% for i=1:nx
%     error(i)=abs((u_analit(i,nt)-u(i,nt))/u(i,nt)) ;  
% end
% 
% figure()
% plot(0.05:0.05:1,error(1:20,1),'k')
% xlim([0.05 1])
% title('Error')

for k=1:nt
    for i=1:nx 
        error(i,k)=abs(u_analit(i,k)-u(i,k));  
    end
end


figure()
contourf(error,100,'linecolor','non');
colormap(jet(256));
colorbar
%caxis([0.1,0.1])
title('error')