%% Homework 6: Question 1
% Write a code to solve the following PDE. If using SOR, I suggest using an 
% over-relaxation factor of 1.5.

% Solve ∇^2(P) = R(x,y) on a unit square domain with zero-gradient BC's.

% To test your code, use the following right-hand side:
% R(x,y) = -4π^2[cos(2πx))+cos(2πy)]

% Note that before you call your Poisson solver, your code should make sure 
% that the integral of the right-hand side is exactly zero, making tiny 
% adjustments as necessary.

% Turn in your code, a contour plot of your solution, the exact solution, 
% and an explanation of why the right-hand side must sum to zero. Also vary 
% the problem size from 20 by 20 points to 100 by 100 and plot, on log-log 
% axes, the observed variation in CPU cost as a function of the number of 
% unknowns.

tic;

N1 = 20;
[p Rtest, x, y] = myPoissonSolver(N1);
toc;

tic;
N2 = 100;
[p2 Rtest2, x2, y2] = myPoissonSolver(N2);
toc;

%% Plot solution first for 20 then for 100

% Plot RHS equation used to test code
figure(1)
contourf(x,y,Rtest')
title('Exact Solution, N = 20');
xlabel('x')
ylabel('y')

% plot solver
figure(2)
contourf(x,y,p');
title('Poisson Solver, N = 20');
xlabel('x')
ylabel('y')

% Plot RHS equation used to test code
figure(3)
contourf(x2,y2,Rtest2')
title('Exact Solution, N=100');
xlabel('x')
ylabel('y')

figure(4)
contourf(x2,y2,p2');
title('Poisson Solver, N = 100');
xlabel('x')
ylabel('y')
%% Create a poisson solver function
function [p, Rtest, x, y] = myPoissonSolver(N)
    % DEFINE DOMAIN
    N = N;     % define number of points in NxN domain
    
    dx = 1/(N-1);
    dy = dx; 
    
    x = 0:dy:1;        % define x with N points
    y = 0:dx:1;        % define y with N points 
    
    % Define right hand side and plug in x and y vectors for Rtest
    R = @(x,y) -4*pi^2*(cos(2*pi*x)+cos(2*pi*y));
    Rtest = zeros(N,N);
    
    % Check that the integral of rhs = 0
    integral = integral2(R,x(1),x(N),y(1),y(N));
    disp("The integral of the right hand side is = "+ integral);
    
    % Redefine R matrix as a vector called RHS
    RHS = zeros(N*N,1);
    
    % Populate Rtest matrix using equation for R and rewrite into RHS
    for i = 1:N
        for j = 1:N
            Rtest(i,j) = R(x(i),y(j));
            RHS((i-1)*N+j) = Rtest(i,j);        % RHS vector 
        end
    end
    
    % SET INTERNAL MATRICCES
    % primary 3 diagonals
    e = ones(N,1);
    e = e.*1/dx^2;
    A_diag = spdiags([e -4*e e], -1:1, N, N);
    
    % off diagonal matrices
    A_off = spdiags([e], 0, N, N);
    
    % Arrange small matrices in large matrix
    A = sparse((N*N), (N*N)); % define sparse matrix, A, which is NxN by NxN
    
    for i = 1:N
        A((i-1)*N+1:(i-1)*N+N,(i-1)*N+1:(i-1)*N+N) = A_diag;        % similar to rearranging the RHS matrix into a vector
    end
    
    for i = 2:N
        A((i-2)*N+1:(i-2)*N+N,(i-1)*N+1:(i-1)*N+N) = A_off;         % top off diagonals 
        A((i-1)*N+1:(i-1)*N+N,(i-2)*N+1:(i-2)*N+N) = A_off;         % bottom off diagonals
    end
    
    % SOLVE FOR P
    p = -A\RHS;
    p = reshape(p, N, N);
end 

%% Would then implement zero boundary conditions... 
% I would rewrite the above matrices to have (N-1)x(N-1) internal small
% matrices and then put those into a large internal matrix and put boundary
% conditions around - which are described in the attached handwritten
% sheet. 
%
% The 100x100 is about 3-5x slower than the 20x20 solver. 
% 
% In both cases the integral of the RHS is very close to 0. If I remember
% correctly in one of the first homeworks we were asked to use Gauss
% Divergence theorem on a poisson eqn and proved that the integral of the
% rhs must be 0. Otherwise we would not have conservation and instead have
% generation inside the control volume. 