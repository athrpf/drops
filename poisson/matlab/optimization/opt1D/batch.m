function batch(index)

% ihcp(tf,L,ntimes,N,q_app,exflux,sigma,nmax)

% tf - final time
% L - length of space intervall
% ntimes - number of temporal steps from t=0 to t=tf
% N - number of meshpoints in space
% q_app - start approximation for heat flux
% exflux - exact heat flux
% sigma - noise
% nmax - max number of iterations

switch index
case 1
    % heat flux without noise
    clear;
    tf= 1.4
    %tf= 7.0
	L= 1
	ntimes= 100
    %ntimes= 500
	N= 10
	q_app= zeros(1,ntimes+1);
	dt= tf/ntimes
    index= 1;
	for i= 0:dt:tf
        exflux(1,index)= q(i);
        index= index+1;
	end
 	sigma= 0
    nmax= 100
    tol= 10^(-6)
	sol= ihcp(tf,L,ntimes,N,q_app,exflux,sigma,nmax,tol);
    % method of nested iterations
%     tol= 10^(-5)
%     sol1= ihcp(tf,L,ntimes,N,q_app,exflux,sigma,nmax,tol);
%     pause;
%     N= 100
%     nmax= 200
%     tol= 10^(-6)
%     q_app= sol1;
%     sol2= ihcp(tf,L,ntimes,N,q_app,exflux,sigma,nmax,tol);
case 2
    % heat flux without noise
    % and different start approximation
    clear;
    tf = 1.4
	L = 1
	ntimes = 100
	N = 10
	q_app = 0.1*ones(1,ntimes+1);
	dt=tf/ntimes;
	index=1;
	for i=0:dt:tf
        exflux(1,index)=q(i);
        index=index+1;
	end
	sigma = 0
    nmax= 10
    tol= 10^(-8)
	sol= ihcp(tf,L,ntimes,N,q_app,exflux,sigma,nmax,tol);
case 3
    % heat flux with noise
    clear;
    tf = 1.4
	L = 1
	ntimes = 100
	N = 10
	q_app = 0.1*ones(1,ntimes+1);
	dt=tf/ntimes;
	index=1;
	for i=0:dt:tf
        exflux(1,index)=q(i);
        index=index+1;
	end
	sigma = 0.1
    nmax= 10
    tol= 10^(-8)
	sol= ihcp(tf,L,ntimes,N,q_app,exflux,sigma,nmax,tol);
case 4
    % heat flux with noise
    % and different number of mesh points
    clear;
    tf = 1.4
	L = 1
	ntimes = 100
	N = 100
	q_app = 0.1*ones(1,ntimes+1);
	dt=tf/ntimes;
	index=1;
	for i=0:dt:tf
        exflux(1,index)=q(i);
        index=index+1;
	end
	sigma = 0.1
    nmax= 10
    tol= 10^(-8)
	sol= ihcp(tf,L,ntimes,N,q_app,exflux,sigma,nmax,tol);
case 5
    % constant heat flux without noise
    clear;
    tf = 1.4
	L = 1
	ntimes = 100
	N = 10
	q_app = 0.1*ones(1,ntimes+1);
	exflux = ones(1,ntimes+1);
    sigma = 0
    nmax= 100
    tol= 10^(-8)
	sol= ihcp(tf,L,ntimes,N,q_app,exflux,sigma,nmax,tol);
case 6
    % constant heat flux with noise
    clear;
    tf = 1.4
	L = 1
	ntimes = 100
	N = 10
	q_app = 0.1*ones(1,ntimes+1);
	exflux = ones(1,ntimes+1);
    sigma = 0.1
    nmax= 10
    tol= 10^(-8)
	sol= ihcp(tf,L,ntimes,N,q_app,exflux,sigma,nmax,tol);
case 7
    % constant heat flux without noise
    % and different length L
    clear;
    tf = 1.4
	L = 0.001
	ntimes = 100
	N = 10
	q_app = 0.99*ones(1,ntimes+1);
	exflux = ones(1,ntimes+1);
    sigma = 0
    nmax= 10
    tol= 10^(-8)
	sol= ihcp(tf,L,ntimes,N,q_app,exflux,sigma,nmax,tol);
case 8
    % constant heat flux without noise
    % and more time steps
    clear;
    tf = 1.4
	L = 1
	ntimes = 140
	N = 100
	q_app = 0*ones(1,ntimes+1);
	dt=tf/ntimes
	index=1;
	for i=0:dt:tf
        exflux(1,index)=q(i);
        index=index+1;
	end
	sigma = 0
    nmax= 10
    tol= 10^(-8)
	sol= ihcp(tf,L,ntimes,N,q_app,exflux,sigma,nmax,tol);
end

return




