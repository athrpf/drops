% claculate and plot temperature over time

% Geometrieparameter
gl= [C_data(1),C_data(2),C_data(3)];

% Anzahl der Intervalle pro Koordinatenrichtung nach
% uniformer Verfeinerung
ni= [C_data(4),C_data(5),C_data(6)];

% Zeitschrittweite
dt= C_data(7);

% Anzahl der Zeitschritte
ndt= C_data(8);

% %for t=1:ndt+1
%     t=200;
%     qh= qcfun(gl,ni,ndt+1,dt,8);
%     plotq(gl,ni,qh,t,1,0)
%     fqc(t)=getframe(gcf);
% %end
%
% return

%--------------- plot temperature over time -------------------

figure(3)

flag= 0;
% flag=0 - Ausgabe ueber y-z-Ebene (x=1)
% flag=1 - Ausgabe fuer y=z=0.5 (x=1)

n=ndt; % Anzahl Zeitschritte Ausgabe

if flag==0
    for t=1:n+1
    %t=ndt+1;
        plot_T(gl,ni,Ti,t,1,flag)
        title('Temperaturverteilung (Querschnitt)');
        fqc(t)=getframe(gcf);
    end
else
    count=1;

    niy=ni(2);
    niz=ni(3);
    npy=niy+1;
    npz=niz+1;

    for t=1:n+1
        % create time-independent matrix
        Tm=reshape(Ti(:,t),npy,npz);

        % output in the middle of the y-z-plane
        Tmp(1,count)=Tm(niy/2+1,niz/2+1);
        count=count+1;
    end
    time=[0:dt:dt*n];
    plot(time,Tmp);
    Tmp(1,ndt+1)
    hold on
    clear Tm;
    clear Tmp;
end

%         hold on
%         plot_T(gl,ni,T,t,1,flag)
%         title('u_h(0.5,0.5,1,t)');
%         hold off

return