%% Задание данных
clc
clear
format long                        % задаём формат на вывод
epsilon = 0.0001;
firstApproach = 0.9;
secondApproach = 0.7;
local = [0,1];

[x1,iter1,time1] = HalfDivision(local,epsilon);
[x2,iter2,time2,k2] = MetodNewton(local, firstApproach, epsilon);
[x3,iter3,time3,k3] = MetodEasyNewton(local, firstApproach, epsilon);
[x4,iter4,time4,k4] = MetodSecant(local, firstApproach, secondApproach, epsilon);

%% Построение графиков
x = 0:0.01:1;
y = targetFunc(x);

tiledlayout(2,2)

%График для метода половинного деления  
nexttile;
x_dot = x1;
[n,m] = size(x_dot);
x1 = x1(m)
iter1
time1
y_dot = zeros(1,m);

plot(x,y);
hold on
scatter(x_dot,y_dot,15,'filled')
grid on
title('Метод половинного деления')
xlabel('Аргумент x')
ylabel('Значение y')

%График для метода Ньютона
nexttile;
x_dot = x2;
[n,m] = size(x_dot);
x2 = x2(m)
iter2
time2
y_dot = zeros(1,m);

plot(x,y, 'LineWidth', 2);
grid on
title('Метод Ньютона')
xlabel('Аргумент x')
ylabel('Значение y')
hold on
scatter(x_dot,y_dot,15,'filled')

for i = 2:iter2+1
    x_kas = x_dot(i):0.00001:x_dot(i-1);
    y_kas = (x_kas - x_dot(i))*k2(i-1);
    plot(x_kas,y_kas)
end

%График для упрощённого метода Ньютона
nexttile;
x_dot = x3;
[n,m] = size(x_dot);
x3 = x3(m)
iter3
time3
y_dot = zeros(1,m);

plot(x,y, 'LineWidth', 2);
grid on
title('Упрощённый метод Ньютона')
xlabel('Аргумент x')
ylabel('Значение y')
hold on
scatter(x_dot,y_dot,15,'filled')

for i = 2:iter3+1
    x_kas = x_dot(i):0.000001:x_dot(i-1);
    y_kas = (x_kas - x_dot(i))*k3(i-1);
    plot(x_kas,y_kas)
end

%График для метода секущих
nexttile;
x_dot = x4;
[n,m] = size(x_dot);
x4 = x4(m)
iter4
time4
y_dot = zeros(1,m);

plot(x,y, 'LineWidth', 2);
grid on
title('Метод секущих')
xlabel('Аргумент x')
ylabel('Значение y')
hold on
scatter(x_dot,y_dot,15,'filled')

for i = 3:iter4+2
    x_kas = x_dot(i):0.000001:x_dot(i-2);
    y_kas = (x_kas - x_dot(i))*k4(i-2);
    plot(x_kas,y_kas)
end

%% Метод половинного деления
function [x,iter,time] = HalfDivision(local, epsilon) 
    a = local(1,1);
    b = local(1,2);
    iter = 0;
    tic
    while (abs(b-a) > 2*epsilon)
        x(iter+1) = (a+b)/2;
        if targetFunc(a)*targetFunc(x(iter+1)) <= 0
            b = x(iter+1);
        else
            a = x(iter+1);
        end
        iter = iter + 1;
    end
    time = toc;
    
    x(iter+1) = (a+b)/2;
end

%% Метод Ньютона
function [x,iter,time,k] = MetodNewton(local, firstApproach, epsilon) 
    a = local(1,1);
    b = local(1,2);
    h = 0.001; %точность шага при дифференцировании
    iter = 0;
    x(iter+1) = firstApproach;
    
    tic
    while(true)
        iter = iter + 1;
        if (x(iter) == b)
            dif = (targetFunc(x(iter))-targetFunc(x(iter)-h))/h;
            k(iter) = dif;
            x0 = x(iter);
            x(iter+1) = x0 - targetFunc(x0)/dif;
        elseif x(iter) == a
            dif = (targetFunc(x+h)-targetFunc(x(iter)))/h;
            k(iter) = dif;
            x0 = x(iter);
            x(iter+1) = x0 - targetFunc(x0)/dif;
        else
            dif = (targetFunc(x(iter)+h)-targetFunc(x(iter)-h))/(2*h);
            k(iter) = dif;
            x0 = x(iter);
            x(iter+1) = x0 - targetFunc(x0)/dif;
        end
        if (abs(x(iter+1)-x0) < epsilon)
            break;
        end
    end
    time = toc;

end

%% Метод Ньютона (упрощённый)
function [x,iter,time,k] = MetodEasyNewton(local, firstApproach, epsilon) 
    a = local(1,1);
    b = local(1,2);
    h = 0.001; %точность шага при дифференцировании
    iter = 0;
    x(iter+1) = firstApproach;
    
    tic
    if x(iter+1) == b
        dif = (targetFunc(x(iter+1))-targetFunc(x(iter+1)-h))/h;
    elseif x(iter+1) == a
        dif = (targetFunc(x(iter+1)+h)-targetFunc(x(iter+1)))/h;
    else
        dif = (targetFunc(x(iter+1)+h)-targetFunc(x(iter+1)-h))/(2*h);
    end
    while(1)
        iter = iter + 1;
        k(iter) = dif;
        x0 = x(iter);
        x(iter+1) = x0 - targetFunc(x0)/dif;
        if (abs(x(iter+1)-x0) < epsilon)
            break;
        end
    end
    time = toc;

end

%% Метод Секущих
function [x,iter,time,k] = MetodSecant(local, firstApproach, secondApproach, epsilon) 
    a = local(1,1);
    b = local(1,2);
    h = 0.001; %точность шага при дифференцировании
    iter = 0;
    x(iter+2) = secondApproach;
    x(iter+1) = firstApproach;
    x0 = x(iter+1);

    tic
    while(1)
        iter = iter + 1;
        k(iter) = (targetFunc(x(iter+1))-targetFunc(x0))/(x(iter+1)-x0);
        help = x(iter+1) - 1/k(iter)*targetFunc(x(iter+1));
        x0 = x(iter+1);
        x(iter+2) = help;
        if (abs(x(iter+2)-x0) < epsilon)
            break;
        end
    end
    time = toc;
end

%% Расчётная функция
function y = targetFunc(x)
    y = 2.^(x-0.1)-1;
    %y = (x-0.3)^3;
end