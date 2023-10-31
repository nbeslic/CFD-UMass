%% Homework 3: Question 5
% You will test three methods of handling advection using finite difference
% methods. The three methods are 1st order upwinding, second-order central 
% differencing, and third-order upwinding (Eqn. 3.12). You will calculate 
% solutions to the one-dimensional linear wave equation on the domain
% from 0 to 1. Use a mesh spacing of 0.01 and a CFL number of 0.15.
% The wave speed is 2.0 and the upwind boundary condition at x=0 is u=0.0
% for all time. The initial condition is a "triangle profile" given. 
% Make a plot of the solution at t = 0.4.

% Domain Setup 
dx = 0.01;      % mesh spacing
x = 0:dx:1;     % domain from 0-1
N = length(x);

% Define CFL/wave speed
CFL = 0.15;
c = 2.0;

% Define time step
dt = CFL*dx/c;                  
t = 0:dt:1;
M = length(t);

% Define upwind boundary condition at x=0 to be u=0
u = zeros(N,M);     % define empty vector u

% Define initial condition of triangle profile
for i=1:N
    if x(i) > 0 && x(i) <= 0.1
        u(i,1) = x(i);
    elseif x(i) > 0.1 && x(i) <= 0.2
        u(i,1) = 0.2 - x(i);
    elseif x(i) > 0.2 && x(i) <=1
        u(i,1) = 0;
    end
end

% Check inital condition by plotting
% figure(1)
% plot(x,u(:,1));
% xlabel('Domain')
% ylabel('u(x)')
% ylim([0,0.2])
% title('Initial Condition')

%% First order upwinding

% create solution vector and give it the initial condition
u1 = u;

% implement scheme with explicit euler timestepping
for n=1:M-1
    for i=2:N
        u1(i,n+1) = -(c*dt/dx)* (u1(i,n) - u1(i-1,n))+ u1(i,n);
    end
end

% plot at t = 0.4
n_val = floor(0.4/dt);

figure(2)
plot(x, u(:,1),'--');         % plot initial condition
hold on
plot(x, u1(:,n_val));           % plot at t = 0.4
xlabel('Domain')
ylabel('u(x,t)')
title('u(x,t) at t=0.4')
ylim([0,0.15])
hold on

%% Second order central differencing

% create solution vector and give it the initial condition
u2 = u;

% implement scheme with explicit euler timestepping
for n=1:M-1
    for i=2:N-1
       % second order differencing scheme
       u2(i,n+1) = u2(i,n) - (c*dt/(dx*2)) * (u2(i+1,n)-u2(i-1,n));
    end
end

plot(x, u2(:,n_val));
hold on

%% Third-order upwinding

% create solution vector and give it the initial condition
u3 = u;

% implement scheme with explicit euler timestepping
for n=1:M-1
    for i=2:N
        if i == 2
            % Implement 1st order upwinding at first point inside domain
            u3(i,n+1) = -(c*dt/dx)* (u1(i,n) - u1(i-1,n))+ u1(i,n);
        elseif i == N
            % Implement 1st order upwinding at last point inside domain
            u3(i,n+1) = -(c*dt/dx)* (u1(i,n) - u1(i-1,n))+ u1(i,n);
        else
            % Third order upwinding scheme
            u3(i,n+1) = u3(i,n) - ((c*dt)/(6*dx))* (2*u3(i+1,n) + 3*u3(i,n) - 6*u3(i-1,n) + u3(i-2,n));
        end
    end
end

plot(x, u3(:,n_val));  
legend('initial condition','first order', 'central difference','3rd order upwinding')
hold off
