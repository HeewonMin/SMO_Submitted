% Yongsu Jung, IDOLAB, Mechanical Engineering, KAIST
% E-mail : yongsu50@kaist.ac.kr

clear all;clc;close all;
% CBS accuracy is set to 10(-4) in Deterministic ATC and 10(-6) in PATC 
%********************************************************************************************

mfilepath = mfilename('fullpath');
pathidx = strfind(mfilepath,filesep);
cdir = mfilepath(1:pathidx(end));
cd(cdir);
addpath([cdir,'dace'])

%********************************Deterministic ATC Setting**********************************************
if ~exist('Det.mat')
    tol = 1.0e-3; % Optimization ÀÇ tolerance
    options=optimoptions('fmincon','GradObj','off','GradConstr','off','TolCon',tol,'TolX',tol,'display','iter','algorithm','sqp');
    ATC_eps = 5*10^(-3);

    c_norm=0;
    x_old(1,:) =[5 5 5 5]; % InitialValue
    NumofLink = 1; % # of Linking variables = # of agrange multipliers
    NumofSystem = 2; % # of subsystems
    v(1,:) = 0*ones(1,NumofLink); % initial Lagrange multiplier
    w(1,:) = 1*ones(1,NumofLink); % initial weight for augmented Lagrangian coordination
    iter=1; % iteration number
    conv=0; % Convergence index
    alpha = ones(1,NumofLink); % step size for parallel solving
    phase = ones(1,NumofLink); % Phase for parallel solving

    %******************* Probabilistic ATC Setting ***************

    
    
    run('SettingInfo_1.m');
    run('SettingInfo_2.m');
    
       

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    while conv==0
          Var_old = []; % variable history.  
          %%%%%%%% Bottom-up ¹æ½Ä %%%%%%%%%%
          SubsystemInfo{1}.iter = iter;
          SubsystemInfo{2}.iter = iter;
          [LinkingOutNew{1},LinkingInNew{1},SharedVarNew{1},VarNew{1},Sample_out{1},Obj_out{1}] = Subsystem(SubsystemInfo{1},v(iter,:),w(iter,:));
          SubsystemInfo{2}.LinkingInVar = LinkingOutNew{1};
          [LinkingOutNew{2},LinkingInNew{2},SharedVarNew{2},VarNew{2},Sample_out{2},Obj_out{2}] = Subsystem(SubsystemInfo{2},v(iter,:),w(iter,:));
          SubsystemInfo{1}.LinkingInPerf = LinkingInNew{2}; 
          for i = 1 : NumofSystem
              SubsystemInfo{i}.Sample_in = Sample_out{i}; 
              SubsystemInfo{i}.Obj_in = Obj_out{i}; 
              SubsystemInfo{i}.Var = VarNew{i};
          end
          Store(iter,:) = [VarNew{1} VarNew{2}];
          
          c = [LinkingOutNew{1} - LinkingInNew{2}]; % Discrepancy(Out-in)
          c_norm(iter+1,:) = norm(c);

          if max([c_norm(iter+1,:)-c_norm(iter,:)]./c_norm(iter+1,:)) < ATC_eps && max(c_norm(iter+1,:))<ATC_eps
                break
          end



          v(iter+1,:) = v(iter,:) +2*w(iter,:).*c;
          x_cur = abs(LinkingOutNew{1});

          w(iter+1,:) = w(iter,:);
          
          c_norm;
          maxnorm(iter) = c_norm(iter+1,1);
%          x_old(iter+1,:) = Var_old;
          iter = iter + 1;

    end

    save('Det.mat','SubsystemInfo', 'x_old', 'options', 'NumofLink', 'NumofSystem');  
    clearvars -except SubsystemInfo x_old options NumofLink NumofSystem
end

load('Det.mat');

%% Probabilistic ATC
%%%%%%%%%%%% Prepare sobol sequence %%%%%%%%%%%%%%%%%%%%%%%%%%%
sn=5*10^5;
qmc=sobolset(1,'Skip',1e3);
qmc=net(qmc,sn);
for i = 1 : NumofSystem
    Info = SubsystemInfo{i};
    for j = 1 : length(Info.Var)-length(Info.Coupledidx)
        SubsystemInfo{i}.Seed{j} = qmc(randperm(length(qmc))); % Sobol sequence 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ATC_eps = 5*10^(-3);
c=0;
c_norm=0;
NumofLink = 1; % # of Linking variables = # of agrange multipliers
NumofSystem = 2; % # of subsystems
v(1,:) = 0*ones(1,NumofLink); % initial Lagrange multiplier
w(1,:) = 1*ones(1,NumofLink); % initial weight for augmented Lagrangian coordination
iter = 1;
conv=0; % Convergence index
alpha = ones(1,NumofLink); % step size 
phase = ones(1,NumofLink); % Phase 

while conv==0

      Var_old = []; % variable history.  
      %%%%%%%% Bottom-up %%%%%%%%%%
      SubsystemInfo{1}.iter = iter;
      SubsystemInfo{2}.iter = iter;
      [LinkingOutNew{1},LinkingInNew{1},SharedVarNew{1},VarNew{1},Sample_out{1},Obj_out{1}] = Subsystem_RBDOwithKDE(SubsystemInfo{1},v(iter,:),w(iter,:));
      SubsystemInfo{2}.LinkingInVar = LinkingOutNew{1};

      [LinkingOutNew{2},LinkingInNew{2},SharedVarNew{2},VarNew{2},Sample_out{2},Obj_out{2}] = Subsystem_RBDOwithKDE(SubsystemInfo{2},v(iter,:),w(iter,:));
      
      SubsystemInfo{1}.LinkingInPerf = LinkingInNew{2}; 
      
      
      for i = 1 : NumofSystem
          SubsystemInfo{i}.Sample_in = Sample_out{i}; 
          SubsystemInfo{i}.Obj_in = Obj_out{i}; 
          SubsystemInfo{i}.Var = VarNew{i};
      end
      
      Store(iter,:) = [VarNew{1} VarNew{2}];
      CDF{iter} = LinkingOutNew{1};
      

      c = [mean(LinkingOutNew{1}(:,1)) - mean(LinkingInNew{2})]; % Discrepancy(Out - In)


      c_norm(iter+1,:) = norm(c)

      if max([c_norm(iter+1,:)-c_norm(iter,:)]./c_norm(iter+1,:)) < ATC_eps && max(c_norm(iter+1,:))<ATC_eps
            break
      end


      v(iter+1,:) = v(iter,:) +2*w(iter,:).*c;
      w(iter+1,:) = w(iter,:);


      c_norm;
      maxnorm(iter) = c_norm(iter+1,1);
      iter = iter + 1;

end







