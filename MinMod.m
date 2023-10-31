%% Homework 6: Question 2
% Solving the 1st order linear wave equation using a minmod flux limiter
% between the first order upwinding scheme and the 2nd order central
% differencing scheme based on a 'triangle' shaped initial condition,
% evaluated at t = 0.4. 

%% Domain Setup 
dx = 0.01;      % mesh spacing
x = 0:dx:1;     % domain from 0-1
I = length(x);

% Define CFL/wave speed
CFL = 0.15;
c = 2.0;

% Define time step
t0 = 0;
tf = 1;
dt = CFL*dx/c;    

t = t0:dt:tf;
N = length(t);

% find value of time t = 0.4
n_val = floor(0.4/dt);

%% Boundary Conditions
% Define upwind boundary condition at x=0 to be u=0
u = zeros(I,N);     % define empty vector u

% Define initial condition of triangle profile
for i=1:I
    if x(i) > 0 && x(i) <= 0.1
        u(i,1) = x(i);
    elseif x(i) > 0.1 && x(i) <= 0.2
        u(i,1) = 0.2 - x(i);
    elseif x(i) > 0.2 && x(i) <=1
        u(i,1) = 0;
    end
end

%% First-order upwinding for comparison

% create solution vector and give it the initial condition
u1 = u;

% implement scheme with explicit euler timestepping
for n=1:N-1
    for i=2:I
        u1(i,n+1) = -(c*dt/dx)* (u1(i,n) - u1(i-1,n))+ u1(i,n);
    end
end

%% Min-Mod Diffusive Flux Limiter 

% Create a solution vector and give it the inital condition
phi = u;

for n=1:N-1
    for i=2:I
        if i == 2
            % first order upwind first interior point
            phi(i,n+1) = -(c*dt/dx)* (phi(i,n) - phi(i-1,n))+ phi(i,n);
        elseif i == I
            % first order upwind last interior point
            phi(i,n+1) = -(c*dt/dx)* (phi(i,n) - phi(i-1,n))+ phi(i,n);
        else
            % min-mod the rest of the interior points
            r = (phi(i,n)-phi(i-1,n))/(phi(i+1,n)-phi(i,n));    % the ratio of two adjacent gradients of the advected scalar
            
            % define min-mod function based on the calculated value of r 
            minmod = max([0, min([r,1])]);

            ud = -(c*dt/dx)* (phi(i,n) - phi(i-1,n))+ phi(i,n);     % upwinding
            cd = phi(i,n) - (c*dt/(dx*2)) * (phi(i+1,n)-phi(i-1,n));    % central differencing
            
            % SCHEME = upwinding - minmod * (upwinding - central differencing)
            phi(i,n+1) = ud - minmod* (ud - cd);
        end
    end
end

%% PLOTTING
% initial condition & 1st order upwinded solution

figure(1)
plot(x, u(:,1),'--');         % plot initial condition
hold on
plot(x, u1(:,n_val));           % plot at t = 0.4
xlabel('Domain')
ylabel('u(x,t)')
title('u(x,t) at t=0.4')
ylim([0,0.15])
hold on 
plot(x, phi(:,n_val));         % plot min-mod solution at t = 0.4
legend('initial condition','1st order upwinding','Min-mod')