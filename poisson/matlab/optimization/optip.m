clear;

% get common data
[C_data,func,niter_max,epsT,M_data,D_data,A_data,S_data]= getParameters;

% Geometrieparameter
gl= [C_data(1),C_data(2),C_data(3)];

% Anzahl der Intervalle pro Koordinatenrichtung nach
% uniformer Verfeinerung
ni= [C_data(4),C_data(5),C_data(6)];

% nix= number of intervals in x-direction
% npx= number of points in x-direction
npx= ni(1)+1;
npy= ni(2)+1;
npz= ni(3)+1;
npyz= npy*npz;
dyd= gl(2)/ni(2);
dzd= gl(3)/ni(3);
dyz= 1e-6*dyd*dzd;

% Zeitschrittweite
dt= C_data(7);
% Anzahl der Zeitschritte
ndt= C_data(8);

npt= ndt+1;
time= [0:dt:ndt*dt];

% create measurement data
[Val,Td]= simdp(C_data,M_data);
DROPS_Val(1,1)= Val;
qcf= qcfun(gl,ni,npt,func);
[Val,Timess]= simip(qcf,C_data,M_data);
DROPS_Val(1,2)= Val;
DROPS_Val(1,3)= 0;

% start approximation
qc= zeros(npyz,ndt+1);
%qc= qcf+1e-3*rand(npyz,ndt+1);

niter= 0;

save OptData func C_data M_data D_data A_data S_data;

% Optimization
while 1    
    niter= niter+1
    % direct problem
    tic
    [Val,Tihat]= simip(qc,C_data,D_data);
    comptime(niter,1)= toc;
    DROPS_Val(niter+1,1)= Val;
    % calculate defect
    dThat= Tihat-Timess;
    norm_err= norm(dThat,'fro')^2*dt*dyz;
    J(niter)= norm_err

    qc_iter(niter,:,:)= qc;
    That_iter(niter,:,:)= Td+Tihat;
    if (norm_err < epsT) | (niter >= niter_max)
        break
    else
        % adjoint problem
        tic
        [Val,lambda]= simipad(fliplr(2*dThat),C_data,A_data);
        comptime(niter,2)= toc;
        DROPS_Val(niter+1,2)= Val;
        Tad= fliplr(lambda);
        % determine search direction 
        if niter==1
            gamma= 0;
            Pn= Tad;
        else        
            gamma= (norm(Tad,'fro')/norm(Tad_old,'fro'))^2;
            Pn= Tad+gamma*Pn;
        end
        % sensitivity problem
        tic
        [Val,Tisens]= simip(Pn,C_data,S_data);
        comptime(niter,3)= toc;
        DROPS_Val(niter+1,3)= Val;
        % date up approximation 
        beta= sum(diag(dThat'*Tisens))/(norm(Tisens,'fro'))^2;
        qc= qc-beta*Pn;  
        Tad_old= Tad;
    end
    save OptResult qc_iter J DROPS_Val;
end
