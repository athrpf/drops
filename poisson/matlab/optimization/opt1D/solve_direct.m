function u = solve_direct(ic,bcx0,bcx1,tf,L)

% Solve the heat equation using backward Euler

% ntimes - number of temporal steps from t=0 to t=tf
ntimes=size(bcx0,2)-1;
% dt - time increment
dt=tf/ntimes;
% N - number of meshpoints in space
N=size(ic,1);
% dx - space increment
dx=L/(N-1);

% time vector
t=linspace(0,tf,ntimes+1);

% initial condition u(:,t=0)
u = ic;

% discretization matrix
e = ones(N,1);
A = spdiags([e -2*e e], -1:1, N, N);

% set boundary conditions
A(1,2) = 2;
A(N,N-1) = 2;

dbc = zeros(N,ntimes+1);
dbc(1,:) = bcx0;
dbc(N,:) = bcx1;

for k = 1 : ntimes
	u(:,k+1) = (eye(N)-dt/dx^2*A)\(u(:,k)+2*dt/dx*dbc(:,k+1));
end
% plot temperature near x=0 as a function of time
%figure(4)
%plot(t,u(1,:),'b')
% add temperature near x=L
%hold on
%plot(t,u(N,:),'g')
%legend('Point 1','Point N')
%title 'temperature vs. time at two points'
%hold off
%pause;

return

