% % Initialisation

tspan = [0:0.1:500];
arr1 = repelem(0.2,2);
arr2 = repelem(0.02,5);
y0 = [arr1 arr2 arr1 arr2];
global mu
global lambda

%%

mu = 0.1501
lambda = 1
[t,y] = ode45(@odefun, tspan, y0);
y(end,:)

% % plotting
% col=[0 0 0;0 0 1;0 1 0;0 1 1;1 0 1;1 0 0;1 1 0;1 0.5 0.5;1 153/255 1;102/255 0 51/255;0 204/255 55/255;153/255 1 0;62/255 62/255 0;150/255 60/255 80/255]
% figure(1); cla;
% hold on
% for i = 1:14
%   plot(t,y(:, i), 'o-', 'Color', col(i,:));
% end
% legend('y(1)','y(2)','y(3)','y(4)','y(5)','y(6)','y(7)','y(8)','y(9)','y(10)','y(11)','y(12)','y(13)','y(14)');

plot(t,y,'o-'),;
legend('y(1)','y(2)','y(3)','y(4)','y(5)','y(6)','y(7)','y(8)','y(9)','y(10)','y(11)','y(12)','y(13)','y(14)');


%% RandomCalculations

% opt = 0;
% iter = 500
% mu1 = rand(1,iter);
% sd = 100;
% for(i = 1:iter)
%     lambda = 1;
%     mu = mu1(i);
%     [t,y] = ode45(@odefun, tspan, y0);
%     if(std(y(end,[3:7,10:14])) < sd)
%         opt = i;
%         sd = std(y(end,[3:7,10:14]))
%     end
% end
% mu1(opt)
% mu = mu1(opt);
% [t,y] = ode45(@odefun, tspan, y0);
% y(end,:)
% sd
%%
function dydt = odefun(t,y)
    global mu
    global lambda
    dydt = zeros(14,1);
    dydt = zeros(14,1);
    dydt(1) = mu*y(8) - lambda*5*y(8)*y(1) + mu*y(2);
    dydt(2) = -lambda*y(2) + mu*y(9) + lambda*5*y(8)*y(1) - mu*y(2);
    dydt(3) = mu*y(10) + mu*y(4);
    dydt(4) = -lambda*y(4) + mu*y(11) - mu*y(4) + 2*mu*y(5);
    dydt(5) = -2*lambda*y(5) + mu*y(12) - 2*mu*y(5) + 3*mu*y(6);
    dydt(6) = -3*lambda*y(6) + mu*y(13) - 3*mu*y(6) + 4*mu*y(7);
    dydt(7) = -4*lambda*y(7) + mu*y(14) - 4*mu*y(7);
    dydt(8) = -mu*y(8) - 5*lambda*y(8)*y(8) + mu*y(9);
    dydt(9) = -mu*y(9) + lambda*y(2) + 5*lambda*y(8)*y(8) - mu*y(9);
    dydt(10) = -mu*y(10) - 4*lambda*y(10) + mu*y(11);
    dydt(11) = -mu*y(11) + lambda*1*y(4) - 3*lambda*y(11) + 4*lambda*y(10) - mu*y(11) + 2*mu*y(12);
    dydt(12) = -mu*y(12) + lambda*2*y(5) - 2*lambda*y(12) + 3*lambda*y(11) - 2*mu*y(12) + 3*mu*y(13);
    dydt(13) = -mu*y(13) + lambda*3*y(6) - 1*lambda*y(13) + 2*lambda*y(12) - 3*mu*y(13) + 4*mu*y(14);
    dydt(14) = -mu*y(14) + lambda*4*y(7) + lambda*y(13) - 4*mu*y(14);
end




