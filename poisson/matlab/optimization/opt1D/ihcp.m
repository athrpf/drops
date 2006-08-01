function [sol]=ihcp(tf,L,ntimes,N,q_app,exflux,sigma,nmax,tol)
%function ihcp(tf,L,ntimes,N,q_app,exflux,sigma,nmax,tol)

% tf - final time
% L - length of space intervall
% ntimes - number of temporal steps from t=0 to t=tf
% N - number of meshpoints in space
% q_app - start approximation for heat flux
% exflux - exact heat flux
% sigma - noise
% nmax - max number of iterations
% tol - tolerance for the defect

% time for plot
t= linspace(0,tf,ntimes+1);
% number of iterations
n= 0;

% time step
dt= tf/ntimes;

% values for boundaries
zer= zeros(1,ntimes+1);
zer_x= zeros(N,1);

% initial temperature
T0= zeros(N,1);

% temperature for exact heat flux
Y= solve_direct(T0,exflux(1,:),zer,tf,L);
Y_x1= Y(N,:);

% standard deviation of measurements
% random variable within -0.5..0.5
omega= rand(1,ntimes+1)-0.5;
% perturbed data
Y_x1= Y_x1+sigma*omega;

% startapproximation
%q_app= exflux-res_q;

% criterion functional
%J= zeros(1,nmax);

while 1

	% loese das direkte Problem
	T= solve_direct(T0,q_app(1,:),zer,tf,L);
    T_x1= T(N,:);

    % plot temperature
    figure(4)
    plot(t,T_x1,'b')
    hold on
    plot(t,Y_x1,'g')
    legend('approximation','exact value')
    title 'temperature vs. time'
    hold off

    % calculate defect
    res= T_x1-Y_x1;
    defect= dt*norm(res)^2;

    J(1,n+1)= defect;
    q_iter(n+1,:,:)= q_app;

    % save data
    save OptResult.mat J q_iter tf L ntimes N exflux sigma tol;

    if (defect >= tol & n<nmax)

        ad_x1= fliplr(2*res);

        % solve adjoint problem
		ad_sol= solve_direct(zer_x,zer,ad_x1(1,:),tf,L);
        lambda= fliplr(ad_sol);

        % determine gradient of functional
		Grad_new= lambda(1,:);
        % correction for the last time step
        %Grad_new(1,end) = Grad_new(1,end-1);
		
        % determine direction of descent with conjugate coefficients
		if n==0
            d= Grad_new;
		else
            gamma= (norm(Grad_new)^2)/(norm(Grad_old)^2);
            d= Grad_new+(gamma*d);
		end
		
        % solve sensitivity problem
		delta_q= d;
		delta_T= solve_direct(zer_x,delta_q(1,:),zer,tf,L);
		delta_Tx1= delta_T(N,:);

        %beta= (norm(Grad_new)^2)/(norm(delta_Tx1)^2);
        beta= dot(res',delta_Tx1')/(norm(delta_Tx1)^2);

        % determine new approximation
        q_app= q_app-beta*d;
    else
        break
    end
    Grad_old= Grad_new;
    n= n+1;

    % plot approximation of heat flux near x=0
 	figure(1)
    plot(t,q_app,'b')

    % plot(t,q_app,'--bs','LineWidth',1,...
    %    'MarkerEdgeColor','k',...
    %    'MarkerFaceColor','w',...
    %    'MarkerSize',2)

	% add exact heat flux
	hold on
	plot(t,exflux,'g')
	legend('approximation','exact value')
	title 'heat flux vs. time'
	hold off

    % plot residual
%     figure(2)
%     res_q= exflux-q_app;
%     plot(t,res_q,'b')
%     title 'residual vs. time'
end

%res_q= exflux-q_app;
%defect_q= norm(res_q)^2;
%display(defect_q);

%display(n);
%display(defect);

figure(3)
semilogy(J, 'b')
title 'cost function J'

sol= q_app;
return










