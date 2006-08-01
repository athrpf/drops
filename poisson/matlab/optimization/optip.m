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

%----- create measurement data - begin

tic
[Val1,Td]= simdp(C_data,M_data);
comptime(1,1)= toc;
DROPS_Iter(:,1)= Val1;

% exakte Messdaten
qcf= qcfun(gl,ni,npt,func);
[Val1,Timess]= simip(qcf,C_data,M_data);

% % Messfehler
% sigma= 0.7;
% omega= randn(size(Timess,1), size(Timess,2));
% Timess= Timess+sigma*omega;

% % bereite die echten Messdaten vor
% R= produceMeasData(100,200,97,197,50);
%
% Rresh= reshape(R,19109,50);
% T0_ls= Rresh(:,1);
%
% Timess=[T0_ls,Rresh]-Td;
% clear R Rresh T0_ls;

% start approximation
qc= 278.652e-3*ones(npyz,ndt+1);

tic
[Val1,Tihat]= simip(qc,C_data,D_data);
comptime(1,2)= toc;
DROPS_Iter(:,2)= Val1;

% Sicherung der Messdaten
Tm= Timess+Td;

%----- create measurement data - end

save OptData func C_data M_data D_data A_data S_data;
save OptResult0 Tm DROPS_Iter comptime

niter= 0;

% Optimization
while 1
    niter= niter+1

    % direct problem
    %Tihat= Tihat-beta*Tisens;

    % calculate defect
    dThat= Tihat-Timess;
    norm_err= norm(dThat,'fro')^2*dt*dyz;
    J= norm_err

    % calculate temperature
    T= Tihat+Td;

    savename= sprintf('OptResult%d',niter);
    save(savename,'qc','T','J','DROPS_Iter','comptime');

    if (norm_err < epsT) | (niter >= niter_max)
        break
    else
        % adjoint problem
        tic
        [Val1,lambda]= simipad(fliplr(2*dThat),C_data,A_data);
        comptime(1,1)= toc;
        DROPS_Iter(:,1)= Val1;
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
        [Val1,Tisens]= simip(Pn,C_data,S_data);
        comptime(1,2)= toc;
        DROPS_Iter(:,2)= Val1;
        % update approximation
        beta= sum(diag(dThat'*Tisens))/(norm(Tisens,'fro'))^2;
        qc= qc-beta*Pn;
        Tad_old= Tad;

        % direct problem
        Tihat= Tihat-beta*Tisens;
    end

end
