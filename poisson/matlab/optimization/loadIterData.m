clear
%cd longIter0_25;
%cd lx3;
disp('load data...');
load OptData;
count= 1;
for iter= 1:20
  filename= sprintf('OptResult%d',iter);
  load(filename);
  %DROPS_Iter_iter(count,:)= DROPS_Iter;
  comptime_iter(count,:,:)= comptime;
  J_iter(count)= J;
  T_iter(count,:,:)= T;
  qc_iter(count,:,:)= qc;
  count= count+1;
end
clear DROPS_Iter J T qc iter filename count comptime;
J(:,:)= J_iter;
%DROPS_Iter= DROPS_Iter_iter;
comptime= comptime_iter;
clear DROPS_Iter_iter J_iter comptime_iter;
disp('complete !');

% % Skalierung des Waermestroms
% qc_iter= 0.023*qc_iter;


% clear
% cd nestedIter;
%
% level= 0;
%
% disp('load data...');
% filename= sprintf('OptData_l%d',level);
% load (filename);
% count= 1;
% for iter= 1:400
%   filename= sprintf('OptResult%d_l%d',iter,level);
%   load(filename);
%   %DROPS_Iter_iter(count,:)= DROPS_Iter;
%   comptime_iter(count,:,:)= comptime;
%   J_iter(count)= J;
%   T_iter(count,:,:)= T;
%   qc_iter(count,:,:)= qc;
%   count= count+1;
% end
% clear DROPS_Iter J T qc iter filename count comptime;
% J(:,:)= J_iter;
% %DROPS_Iter= DROPS_Iter_iter;
% comptime= comptime_iter;
% clear DROPS_Iter_iter J_iter comptime_iter;
% disp('complete !');
%
% % % Skalierung des Waermestroms
% % qc_iter= 0.023*qc_iter;