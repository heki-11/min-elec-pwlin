%%---------------------- ME 558 Project -----------------------------%%
% Real-time pricing!
% Converts nonlinear (but piecewise linear) objective into objective by
% using slack variables
% Assumes quadratic relationship between price and load (i.e., price =
% (load)^2)

clear all
close all
clc

%% Define parameters

T = 1/2*60; % Total time (min)
m = 3; % Number of machines
N0 = 5; % Predefined quota
a = 3;
b = 3.6052;

% Processing time (min)
p1 = 3;
p2 = 4;
p3 = 3;

q1_pk_pr = 7;
q2_pk_pr = 7;
q3_pk_pr = 9;
q = [q1_pk_pr q2_pk_pr q3_pk_pr];

tot_len = 3*m*T + 4*m - 1;
tot_len_mod = tot_len + T + 1 + (T+1)+5;

%% Equality constraints

% Constraint 2

Aeq2 = zeros(1,tot_len_mod);
ind = 1;

i = 1;
    for j = 0:p1-1
        k = (i-1)*(T+1) + (j+1);
        Aeq2(ind,k) = 1;
        beq2(ind,1) = 0;
        ind = ind + 1;
    end
i = 2;
    for j = 0:p2-1
        k = (i-1)*(T+1) + (j+1);
        Aeq2(ind,k) = 1;
        beq2(ind,1) = 0;
        ind = ind + 1;
    end
i = 3;
    for j = 0:p3-1
        k = (i-1)*(T+1) + (j+1);
        Aeq2(ind,k) = 1;
        beq2(ind,1) = 0;
        ind = ind + 1;
    end

   

% Constraint 3

Aeq3 = zeros(1,tot_len_mod);
ind = 1;

i = 1;
    for j = p1:T
        k = (i-1)*(T+1) + (j+1);
        k1 = (i-1)*(T+1) + (0+1) + 2*m*(T+1);
        k2 = (i-1)*(T+1) + (j-p1+1+1) + 2*m*(T+1);
        Aeq3(ind,k) = 1;
        Aeq3(ind,k1:k2) = -1;
        beq3(ind,1) = 0;
        ind = ind + 1;
    end
i = 2;
    for j = p2:T
        k = (i-1)*(T+1) + (j+1);
        k1 = (i-1)*(T+1) + (0+1) + 2*m*(T+1);
        k2 = (i-1)*(T+1) + (j-p2+1+1) + 2*m*(T+1);
        Aeq3(ind,k) = 1;
        Aeq3(ind,k1:k2) = -1;
        beq3(ind,1) = 0;
        ind = ind + 1;
    end
i = 3;
    for j = p3:T
        k = (i-1)*(T+1) + (j+1);
        k1 = (i-1)*(T+1) + (0+1) + 2*m*(T+1);
        k2 = (i-1)*(T+1) + (j-p3+1+1) + 2*m*(T+1);
        Aeq3(ind,k) = 1;
        Aeq3(ind,k1:k2) = -1;
        beq3(ind,1) = 0;
        ind = ind + 1;
    end

% Constraint 9

Aeq9 = zeros(1,tot_len_mod);
ind = 1;

i = 1;
    for j = 1:p1-1
        k = (i-1)*(T+1) + (j+1) + m*(T+1);
        k1 = (i-1)*(T+1) + (0+1) + 2*m*(T+1);
        k2 = (i-1)*(T+1) + (j+1) + 2*m*(T+1);
        Aeq9(ind,k) = 1;
        Aeq9(ind,k1:k2) = -1;
        beq9(ind,1) = 0;
        ind = ind + 1;
    end
i = 2;
    for j = 1:p2-1
        k = (i-1)*(T+1) + (j+1) + m*(T+1);
        k1 = (i-1)*(T+1) + (0+1) + 2*m*(T+1);
        k2 = (i-1)*(T+1) + (j+1) + 2*m*(T+1);
        Aeq9(ind,k) = 1;
        Aeq9(ind,k1:k2) = -1;
        beq9(ind,1) = 0;
        ind = ind + 1;
    end
i = 3;
    for j = 1:p3-1
        k = (i-1)*(T+1) + (j+1) + m*(T+1);
        k1 = (i-1)*(T+1) + (0+1) + 2*m*(T+1);
        k2 = (i-1)*(T+1) + (j+1) + 2*m*(T+1);
        Aeq9(ind,k) = 1;
        Aeq9(ind,k1:k2) = -1;
        beq9(ind,1) = 0;
        ind = ind + 1;
    end


% Constraint 10

ind = 1;
Aeq10 = zeros(1,tot_len_mod);

i = 1;
    for j = p1:T
        k = (i-1)*(T+1) + (j+1) + m*(T+1);
        k1 = (i-1)*(T+1) + (j-p1+1+1) + 2*m*(T+1);
        k2 = (i-1)*(T+1) + (j+1) + 2*m*(T+1);
        Aeq10(ind,k) = 1;
        Aeq10(ind,k1:k2) = -1;
        beq10(ind,1) = 0;
        ind = ind + 1;
    end
i = 2;
    for j = p2:T
        k = (i-1)*(T+1) + (j+1) + m*(T+1);
        k1 = (i-1)*(T+1) + (j-p2+1+1) + 2*m*(T+1);
        k2 = (i-1)*(T+1) + (j+1) + 2*m*(T+1);
        Aeq10(ind,k) = 1;
        Aeq10(ind,k1:k2) = -1;
        beq10(ind,1) = 0;
        ind = ind + 1;
    end
i = 3;
    for j = p3:T
        k = (i-1)*(T+1) + (j+1) + m*(T+1);
        k1 = (i-1)*(T+1) + (j-p3+1+1) + 2*m*(T+1);
        k2 = (i-1)*(T+1) + (j+1) + 2*m*(T+1);
        Aeq10(ind,k) = 1;
        Aeq10(ind,k1:k2) = -1;
        beq10(ind,1) = 0;
        ind = ind + 1;
    end



%% Inequality constraints

% Constraint 4

ind = 1;
Aeq4 = zeros(1,tot_len_mod);

i = 1;
    for j = 0:T-p3-p2
        k = (i-1)*(T+1) + (j+1);
        k1 = (i+1-1)*(T+1) + (j+1);
        k2 = (i+1-1)*(T+1) + (j+1) + m*(T+1);
        Aeq4(ind,k1) = 1;
        Aeq4(ind,k) = -1;
        Aeq4(ind,k2) = 1;
        beq4(ind,1) = 0;
        ind = ind + 1;
    end

i = 2;
    for j = 0:T-p3
        k = (i-1)*(T+1) + (j+1);
        k1 = (i+1-1)*(T+1) + (j+1);
        k2 = (i+1-1)*(T+1) + (j+1) + m*(T+1);
        Aeq4(ind,k1) = 1;
        Aeq4(ind,k) = -1;
        Aeq4(ind,k2) = 1;
        beq4(ind,1) = 0;
        ind = ind + 1;
    end

% Constraint 5

ind = 1;
Aeq5 = zeros(1,tot_len_mod);

for i = 1:m-1
    for j = 0:T
        k = (i-1)*(T+1) + (j+1);
        k1 = (i+1-1)*(T+1) + (j+1);
        k2 = i + 3*m*(T+1);
        Aeq5(ind,k1) = -1;
        Aeq5(ind,k) = 1;
        Aeq5(ind,k2) = -1;
        beq5(ind,1) = 0;
        ind = ind + 1;
    end
end 

% Constraint 6

i = m;
j = T;
ind = 1;
k = (i-1)*(T+1) + (j+1);
Aeq6 = zeros(1,tot_len_mod);

Aeq6(ind,k) = -1;
beq6(ind,1) = -N0;

% Constraint 8

ind = 1;
Aeq8_1 = zeros(1,tot_len_mod);
Aeq8_2 = zeros(1,tot_len_mod);
Aeq8_3 = zeros(1,tot_len_mod);
Aeq8_4 = zeros(1,tot_len_mod);

for i = 1:m-1
    for j = 0:T
        k = (i-1)*(T+1) + (j+1) + m*(T+1);
        Aeq8_1(ind,k) = 1;
        beq8_1(ind,1) = 1;
        ind = ind + 1;
    end
end 

ind = 1;

for i = 1:m-1
    for j = 0:T
        k = (i-1)*(T+1) + (j+1) + m*(T+1);
        Aeq8_2(ind,k) = -1;
        beq8_2(ind,1) = 0;
        ind = ind + 1;
    end
end 

ind = 1;

for i = 1:m-1
    for j = 0:T
        k = (i-1)*(T+1) + (j+1) + 2*m*(T+1);
        Aeq8_3(ind,k) = 1;
        beq8_3(ind,1) = 1;
        ind = ind + 1;
    end
end 

ind = 1;

for i = 1:m-1
    for j = 0:T
        k = (i-1)*(T+1) + (j+1) + 2*m*(T+1);
        Aeq8_4(ind,k) = -1;
        beq8_4(ind,1) = 0;
        ind = ind + 1;
    end
end 

% Constraint 11

Aeq11 = zeros(1,tot_len_mod);
ind = 1;

i = 1;
    for j = 1:T-p1+1
        k = (i-1)*(T+1) + (j+1) + 2*m*(T+1);
        k1 = (i-1)*(T+1) + (j+1) + m*(T+1);
        k2 = (i-1)*(T+1) + (j+p1-1+1) + m*(T+1);
        Aeq11(ind,k) = p1;
        Aeq11(ind,k1:k2) = -1;
        beq11(ind,1) = 0;
        ind = ind + 1;
    end
i = 2;
    for j = 1:T-p2+1
        k = (i-1)*(T+1) + (j+1) + 2*m*(T+1);
        k1 = (i-1)*(T+1) + (j+1) + m*(T+1);
        k2 = (i-1)*(T+1) + (j+p2-1+1) + m*(T+1);
        Aeq11(ind,k) = p2;
        Aeq11(ind,k1:k2) = -1;
        beq11(ind,1) = 0;
        ind = ind + 1;
    end
i = 3;
    for j = 1:T-p3+1
        k = (i-1)*(T+1) + (j+1) + 2*m*(T+1);
        k1 = (i-1)*(T+1) + (j+1) + m*(T+1);
        k2 = (i-1)*(T+1) + (j+p3-1+1) + m*(T+1);
        Aeq11(ind,k) = p3;
        Aeq11(ind,k1:k2) = -1;
        beq11(ind,1) = 0;
        ind = ind + 1;
    end

%% Constraint on slack

qpos = [0;7;9;14;16;23];

for i = 0:T
    comm(i+1) = -40/T^2*(i - T/2)^2 + 13;
end

for i = 0:T
%         y1(i+1) = exp(-0.4*(i-18.75)^2+200 + a*qpos(1) + b)*qpos(1);
%         y2(i+1) = exp(-0.4*(i-18.75)^2+200 + a*qpos(2) + b)*qpos(2);
%         y3(i+1) = exp(-0.4*(i-18.75)^2+200 + a*qpos(3) + b)*qpos(3);
%         y4(i+1) = exp(-0.4*(i-18.75)^2+200 + a*qpos(4) + b)*qpos(4);
%         y5(i+1) = exp(-0.4*(i-18.75)^2+200 + a*qpos(5) + b)*qpos(5);
%         y6(i+1) = exp(-0.4*(i-18.75)^2+200 + a*qpos(6) + b)*qpos(6);
        y1(i+1) = qpos(1)^2 + comm(i+1);
        y2(i+1) = qpos(2)^2 + comm(i+1);
        y3(i+1) = qpos(3)^2 + comm(i+1);
        y4(i+1) = qpos(4)^2 + comm(i+1);
        y5(i+1) = qpos(5)^2 + comm(i+1);
        y6(i+1) = qpos(6)^2 + comm(i+1);
end

ind = 1;
Aeqslack = zeros(1,tot_len_mod);
for j = 0:T
    k = tot_len + j + 1;
    Aeqslack(ind,k) = 1;
    Aeqslack(ind,tot_len + (T + 1) + j + 1) = -(y2(j+1)-y1(j+1))/(qpos(2)-qpos(1));
    beqslack(ind,1) = y1(j+1);
    
    Aeqslack(ind+1,k) = 1;
    Aeqslack(ind+1,tot_len + (T + 1) + j + 1) = -(y3(j+1)-y2(j+1))/(qpos(3)-qpos(2));
    beqslack(ind+1,1) = y2(j+1);
    
    Aeqslack(ind+2,k) = 1;
    Aeqslack(ind+2,tot_len + (T + 1) + j + 1) = -(y4(j+1)-y3(j+1))/(qpos(4)-qpos(3));
    beqslack(ind+2,1) = y3(j+1);
    
    Aeqslack(ind+3,k) = 1;
    Aeqslack(ind+3,tot_len + (T + 1) + j + 1) = -(y5(j+1)-y4(j+1))/(qpos(5)-qpos(4));
    beqslack(ind+3,1) = y4(j+1);
    
    Aeqslack(ind+4,k) = 1;
    Aeqslack(ind+4,tot_len + (T + 1) + j + 1) = -(y6(j+1)-y5(j+1))/(qpos(6)-qpos(5));
    beqslack(ind+4,1) = y5(j+1);
    
    ind = ind + 5;
    
end

ind = 1;
Aeqslackrel = zeros(1,tot_len_mod);
for j = 0:T
    k = (1-1)*(T+1) + (j+1) + m*(T+1);
    k2 = (2-1)*(T+1) + (j+1) + m*(T+1);
    k3 = (3-1)*(T+1) + (j+1) + m*(T+1);
    k4 = tot_len + T + 1 + j + 1;
    Aeqslackrel(ind,k) = q(1);
    Aeqslackrel(ind,k2) = q(2);
    Aeqslackrel(ind,k3) = q(3);
    Aeqslackrel(ind,k4) = -1;
    beqslackrel(ind,1) = 0;
    ind = ind + 1;
    
end
    



%% Concatenation of constraints

Aeq = [Aeq2;Aeq3;Aeq9;Aeq10;Aeqslackrel];
beq = [beq2;beq3;beq9;beq10;beqslackrel];

Aineq = [Aeq4;Aeq5;Aeq6;Aeq11;-Aeqslack];
bineq = [beq4;beq5;beq6;beq11;-beqslack];

    
%% Lower and upper bound (binary enforcement & positivity)

lb = zeros(tot_len_mod,1);
ub = [(N0+2)*ones(m*(T+1),1);1*ones(2*m*(T+1),1);N0*ones(m-1,1);inf*ones(T+1,1);23*ones((T+1)+5,1)];

%% Objective function


f = zeros(1,tot_len + T + 1 + (T+1)+5);
for j = 0:T
    k = tot_len + j + 1;
    f(1,k) = 1;
end
    
% f(1,(m-1)*(T+1)+1:m*(T+1)) = q4;


%% Optimization

options = optimoptions(@intlinprog,'Display','iter',...
    'Heuristics','advanced');


% sol = intlinprog(f,1:3*m*(T+1),Aineq,bineq,Aeq,beq,lb,ub,[],options);

% Aineq = [Aineq;Aeq;-Aeq];
% bineq = [bineq;beq;-beq];

Aineq = sparse(Aineq);
Aeq = sparse(Aeq);
bineq = sparse(bineq);
beq = sparse(beq);

sol = intlinprog(100000*f',union([1:3*m*(T+1)],[tot_len+T+1+1:tot_len+T+1+(T+1)]),...
    Aineq,bineq,Aeq,beq,lb,ub,[],options);
sol_bl = intlinprog(zeros(size(f')),union([1:3*m*(T+1)],[tot_len+T+1+1:tot_len+T+1+(T+1)]),...
    Aineq,bineq,Aeq,beq,lb,ub,[],options);

N = sol(1:m*(T+1));
X = sol(m*(T+1)+1:2*m*(T+1));
Y = sol(2*m*(T+1)+1:3*m*(T+1));
B = sol(3*m*(T+1)+1:3*m*(T+1)+m-1);

% total_cost = f*sol/1000;

N_bl = sol_bl(1:m*(T+1));
X_bl = sol_bl(m*(T+1)+1:2*m*(T+1));
Y_bl = sol_bl(2*m*(T+1)+1:3*m*(T+1));
B_bl = sol_bl(3*m*(T+1)+1:3*m*(T+1)+m-1);

% total_cost_bl = f*sol_bl/1000;



%% Plot results of proposed

time = 0:1:T;
Q = [];

for i = 0:T
    Q = [Q q];
end

for j = 0:T
    ind1 = (1-1)*(T+1) + (j+1) + m*(T+1);
    ind2 = (2-1)*(T+1) + (j+1) + m*(T+1);
    ind3 = (3-1)*(T+1) + (j+1) + m*(T+1);
    price(j+1) = ((q*[sol(ind1);sol(ind2);sol(ind3)])^2 + comm(j+1)).*...
        (q*[sol(ind1);sol(ind2);sol(ind3)]);
    price_bl(j+1) = ((q*[sol_bl(ind1);sol_bl(ind2);sol_bl(ind3)])^2 + comm(j+1)).*...
        (q*[sol_bl(ind1);sol_bl(ind2);sol_bl(ind3)]);
end

total_cost = sum(price)
total_cost_bl = sum(price_bl)

figure;
subplot(4,1,1)
plot(time,price,'k','linewidth',2)
xlabel('Time [min]','fontweight','bold','fontsize',12)
ylabel('Price','fontweight','bold','fontsize',12)
title('Proposed','fontweight','bold','fontsize',12)
% title(['Proposed: N0 = ' num2str(N0) ', Cost = ' num2str(total_cost) ' cents'],'fontweight','bold','fontsize',12)


subplot(4,1,2)
plot(time,X(1:T+1),'bo-','linewidth',2)
hold on
plot(time,X(T+2:2*T+2),'rx-','linewidth',2)
plot(time,X(2*T+3:3*T+3),'g*-')
grid on
xlabel('Time [min]','fontweight','bold','fontsize',12)
ylabel('X','fontweight','bold','fontsize',12)
lg = legend({'Machine 1','Machine 2','Machine 3'},'location','northwest',...
    'orientation','horizontal');
lg.FontSize = 10;

subplot(4,1,3)
plot(time,Y(1:T+1),'bo-','linewidth',2)
hold on
plot(time,Y(T+2:2*T+2),'rx-','linewidth',2)
plot(time,Y(2*T+3:3*T+3),'g*-')
grid on
xlabel('Time [min]','fontweight','bold','fontsize',12)
ylabel('Y','fontweight','bold','fontsize',12)


subplot(4,1,4)
plot(time,N(1:T+1),'bo-','linewidth',2)
hold on
plot(time,N(T+2:2*T+2),'rx-','linewidth',2)
plot(time,N(2*T+3:3*T+3),'g*-')
grid on
ylabel('N','fontweight','bold','fontsize',12)
xlabel('Time [min]','fontweight','bold','fontsize',12)


%% Plot results of baseline (just feasible)


figure;
subplot(4,1,1)
plot(time,price_bl,'k','linewidth',2)
xlabel('Time [min]','fontweight','bold','fontsize',12)
ylabel('Price','fontweight','bold','fontsize',12)
% title(['Proposed: N0 = ' num2str(N0) ', Cost = ' num2str(total_cost) ' cents'],'fontweight','bold','fontsize',12)
% lg = legend({'Machine 1','Machine 2','Machine 3'},'location','northwest');
% lg.FontSize = 10;
title('Baseline','fontweight','bold','fontsize',12)


subplot(4,1,2)
plot(time,X_bl(1:T+1),'bo-','linewidth',2)
hold on
plot(time,X_bl(T+2:2*T+2),'rx-','linewidth',2)
plot(time,X_bl(2*T+3:3*T+3),'g*-')
grid on
xlabel('Time [min]','fontweight','bold','fontsize',12)
ylabel('X','fontweight','bold','fontsize',12)
lg = legend({'Machine 1','Machine 2','Machine 3'},'location','northwest',...
    'orientation','horizontal');
lg.FontSize = 10;

subplot(4,1,3)
plot(time,Y_bl(1:T+1),'bo-','linewidth',2)
hold on
plot(time,Y_bl(T+2:2*T+2),'rx-','linewidth',2)
plot(time,Y_bl(2*T+3:3*T+3),'g*-')
grid on
xlabel('Time [min]','fontweight','bold','fontsize',12)
ylabel('Y','fontweight','bold','fontsize',12)


subplot(4,1,4)
plot(time,N_bl(1:T+1),'bo-','linewidth',2)
hold on
plot(time,N_bl(T+2:2*T+2),'rx-','linewidth',2)
plot(time,N_bl(2*T+3:3*T+3),'g*-')
grid on
ylabel('N','fontweight','bold','fontsize',12)
xlabel('Time [min]','fontweight','bold','fontsize',12)