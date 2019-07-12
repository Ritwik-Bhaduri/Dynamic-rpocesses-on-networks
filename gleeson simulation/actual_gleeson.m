tspan = [0:0.5:50];

y0 = [0.232 0.176 0.008 0.032 0.032 0.04 0.008 0.248 0.144 0.012 0.016 0.032 0.016 0.004 0 0 ]; 
%random initialisation of each vertex to either infected or succeptible or infected

[t,y] = ode45(@odefun, tspan, y0);

col=[0.5 1 0.5;0 1 1;1 0 1;1 0 0;0 0 0;0 0 1;1 .8 0;0.5 0.5 0.5;1 153/255 1;102/255 0 51/255;0 204/255 55/255;153/255 1 0;62/255 62/255 0;150/255 60/255 80/255];
figure(1); cla;
hold on
for i = 1:14
  plot(t,y(:, i),'-', 'Color', col(i,:),'LineWidth',3);
end
legend('s 1,0','s 1,1','s 4,0','s 4,1','s 4,2','s 4,3','s 4,4','i 1,0','i 1,1','i 4 ,0','i 4,1','i 4,2','i 4,3','i 4,4');

s=0;
for i = 1:14
    s = s + y(:,i);
end
plot(t,s,'-', 'Color', [0 0 0],'LineWidth',3);

function dydt = odefun(t,y)
dydt = zeros(14,1);
mu = 1;
lambda = 5;
bs = lambda * (3 * y(4) + 4 * y(5) + 3 * y(6)) / (4 * y(1) + 4 * y(3) + 3 * y(4) + 2 * y(5) + 1 * y(6));
bi = lambda * (4 * y(2) + y(4) + 4 * y(5) + 9 * y(6) + 16 * y(7)) / (4 * y(2) + y(4) + 2 * y(5) + 3 * y(6) + 4 * y(7));
dydt = zeros(14,1);
dydt(1) = mu*y(8) - bs*y(1) + mu*y(2);
dydt(2) = -lambda*y(2) + mu*y(9) + bs*y(1)- mu*y(2);
dydt(3) = mu*y(10) - 4*bs*y(3) + mu*y(4);
dydt(4) =   -lambda*y(4) + mu*y(11) - 3*bs*y(4) + 4*bs*y(3) - 1*mu*y(4) + 2*mu*y(5);
dydt(5) = -2*lambda*y(5) + mu*y(12) - 2*bs*y(5) + 3*bs*y(4) - 2*mu*y(5) + 3*mu*y(6);
dydt(6) = -3*lambda*y(6) + mu*y(13) - 1*bs*y(6) + 2*bs*y(5) - 3*mu*y(6) + 4*mu*y(7);
dydt(7) = -4*lambda*y(7) + mu*y(14) + 1*bs*y(6) - 4*mu*y(7);
dydt(8) = -mu*y(8) -bi*y(8) + mu*y(9);
dydt(9) = -mu*y(9) + lambda*y(2) + bi*y(8) - mu*y(9);
dydt(10) = -mu*y(10) - 4*bi*y(10) + mu*y(11);
dydt(11) = -mu*y(11) + lambda*1*y(4) - 3*bi*y(11) + 4*bi*y(10) - mu*y(11)+2*mu*y(12);
dydt(12) = -mu*y(12) + lambda*2*y(5) - 2*bi*y(12) + 3*bi*y(11) - 2*mu*y(12)+3*mu*y(13);
dydt(13) = -mu*y(13) + lambda*3*y(6) - 1*bi*y(13) + 2*bi*y(12) - 3*mu*y(13)+4*mu*y(14);
dydt(14) = -mu*y(14) + lambda*4*y(7) + 1*bi*y(13) - 4*mu*y(14);
dydt = [dydt; bs; bi];
end