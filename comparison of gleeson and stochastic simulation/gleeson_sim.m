%tspan = linspace(0, (2* iter * 3 * 10^(-4)),iter*2);

%tspan = linspace(0, max(t), iter/5000);    
% t used here comes from output of gillespie simulation. For independent use the above command

tspan = linspace(0, 6, 6 * 10^4);
k_max = max(sum(G))
sus_0 = zeros(k_max+1,k_max+1);
inf_0 = zeros(k_max+1,k_max+1);


for(i = 1:n)
    degree_i = sum(G(i,:));
    infected_neighbors_of_i = sum(G(i,:) * state.');
    if(state(i) == 0)   
        sus_0(degree_i+1, infected_neighbors_of_i+1) = sus_0(degree_i+1, infected_neighbors_of_i+1) + 1;
    else
        inf_0(degree_i+1, infected_neighbors_of_i+1) = inf_0(degree_i+1, infected_neighbors_of_i+1) + 1;
    end
end

y_0 = [reshape(sus_0.', (k_max+1)^2, 1); reshape(inf_0.', (k_max+1)^2, 1)];
% what this line does is join the susceptible and infected intp one array
% so that it can be put in ode45. Here, the matrices are converted by row.  
y_0 = y_0 / sum(y_0);
% y_0 dentes the the proportions

[t, y] = ode45(@odefun, tspan, y_0);


temp = size(y);
temp = temp(1);
inf = zeros(temp,1);
for (i = ((k_max+1)^2 + 1) : 2*(k_max+1)^2)
  inf = inf + y(:,i);
end

hold on;
lh.Color=[0.5,0.1,0.1,1];
plot(t,inf,'color',lh.Color);



function  dydt = odefun(t, y)
    k_max =20;
    dydt = zeros(2 * (k_max + 1)^2, 1);
    mu =0;
    lambda = 1;
    
    y_s = y(1:(k_max+1)^2, 1);
    y_i = y([((k_max+1)^2 + 1) : (2 * (k_max+1)^2)], 1);
    s = reshape(y_s, k_max + 1, k_max + 1).';
    i = reshape(y_i.', k_max + 1, k_max + 1).';
    
    P = sum(s,2) + sum(i,2);
    
    numerator_beta_s = 0;
    for(k = 0:k_max)
        temp = 0;
        for(m = 0:k)
            temp = temp + (k-m) * lambda * m * s(k+1,m+1);
        end
        numerator_beta_s = numerator_beta_s + temp * P(k+1);
    end
    denominator_beta_s = 0;
    for(k = 0:k_max)
        temp = 0;
        for(m = 0:k)
            temp = temp + (k-m) * s(k+1,m+1);
        end
        denominator_beta_s = denominator_beta_s + temp * P(k+1);
    end
    beta_s = numerator_beta_s / denominator_beta_s;
    
    numerator_beta_i = 0;
    for(k = 0:k_max)
        temp = 0;
        for(m = 0:k)
            temp = temp + m * lambda * m * s(k+1,m+1);
        end
        numerator_beta_i = numerator_beta_i + temp * P(k+1);
    end
    denominator_beta_i = 0;
    for(k = 0:k_max)
        temp = 0;
        for(m = 0:k)
            temp = temp + m * s(k+1,m+1);
        end
        denominator_beta_i = denominator_beta_i + temp * P(k+1);
    end
    beta_i = numerator_beta_i / denominator_beta_i;
    
    gamma_s = mu;
    gamma_i = mu;
    
    dsdt = zeros((k_max + 1), (k_max + 1));
    didt = zeros((k_max + 1), (k_max + 1));
    
    for(k = 0:k_max)
        for(m = 1:(k-1))
            dsdt(k+1,m+1) = -lambda*m*s(k+1,m+1) + mu*i(k+1,m+1) - beta_s*(k-m)*s(k+1,m+1) + beta_s*(k-m+1)*s(k+1,m) - gamma_s*m*s(k+1,m+1) + gamma_s*(m+1)*s(k+1,m+2);
            didt(k+1,m+1) = -mu*i(k+1,m+1) + lambda*m*s(k+1,m+1) - beta_i*(k-m)*i(k+1,m+1) + beta_i*(k-m+1)*i(k+1,m) - gamma_i*m*i(k+1,m+1) + gamma_i*(m+1)*i(k+1,m+2);
        end
    end
    
    for(k = 0:k_max)    % when m = 0
        dsdt(k+1,1) = mu*i(k+1,1) - beta_s*(k)*s(k+1,1) + gamma_s*1*s(k+1,2);
        didt(k+1,1) = -mu*i(k+1,1) - beta_i*(k)*i(k+1,1) + gamma_i*1*i(k+1,2);
    end
    for(k = 1:k_max)    % when m = k
        dsdt(k+1,k+1) = -lambda*k*s(k+1,k+1) + mu*i(k+1,k+1) + beta_s*1*s(k+1,k) - gamma_s*k*s(k+1,k+1);
        didt(k+1,k+1) = -mu*i(k+1,k+1) + lambda*k*s(k+1,k+1) + beta_i*1*i(k+1,k) - gamma_i*k*i(k+1,k+1);
    end
    
    
    dydt = [reshape(dsdt.', (k_max+1) * (k_max+1), 1); reshape(didt.', (k_max+1) * (k_max+1), 1)];
    
end