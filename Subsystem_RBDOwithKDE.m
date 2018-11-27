%%% General Subsystem in ATC %%%%
% Yongsu Jung, IDOLAB, Mechanical Engineering, KAIST
% E-mail : yongsu50@kaist.ac.kr
%% Input 정보
% Var : 이전 iteration 의 design variable.
% Linking : 이전 iteration 의 Linking variable.
% ObjFile : 목적함수 function.
% ConstFile : 제한조건 function.
% FEAFile : 제한조건 및 Linking variable 계산함수.
% v : Lagrange multiplier - 모든 Lagrange multiplier 를 포함하므로 idx 로 적절한 index.
% w : Weight.
% Sample_in : 축적된 샘플링 정보.
% Obj_in : 축적된 샘플링 정보.(제한조건)
% Linking_in ; 축적된 샘플링 정보.(Coupling Variable)
% Bound : Design variable 의 Upper and lower bound - Two row matrix.
% Options : fmincon 에 들어가는 option.
% idx : 현재 subsystem 에 사용되는 Lagrange multiplier 에 대한 index.
%% Output 정보
% LinkingNew : 이번 iteration 이 끝난 뒤 linkingvariable 정보 - Distribution 으로 전달해야 한다.
% VarNew : 이번 iteration 이 끝난 뒤 design variable 정보.
% Sample_out ; 이번 iteration 이 끝난 뒤 Sampling 정보.
% LinkingIndex : Linking variable 중이 Performance function 의 몇 번째부터 위치하는가?
% SharedIndex : Linking variable 중에서 variable 중 어떤 variable 들이 linking
% variable 인가?


function [LinkingOutNew,LinkingInNew,SharedVarNew,VarNew,Sample_out,Obj_out] = Subsystem_RBDOwithKDE(Info,v,w)

Globaliter = Info.iter;

Localidx = Info.Localidx; % Var 중에 어떤 것이 Local variable 인지
Coupledidx = Info.Coupledidx; % Var 중에 어떤 것이 Coupled variable 인지
Sharedidx = Info.Sharedidx; % Var 중에 어떤 것이 Shared variable 인지

SharedVarStd = Info.SharedVarStd; % PATC 에서 쓰일 Shared variable 의 std 값.
SharedVarType = Info.SharedVarType; % PATC 에서 쓰일 Shared variable 의 분포

LocalVarStd = Info.LocalVarStd; % PATC 에서 쓰일 local variable 의 std 값.
LocalVarType = Info.LocalVarType; % PATC 에서 쓰일 local variable 의 분포

Var = Info.Var;

LinkingPerfidx = Info.LinkingPerfidx; % FEA analysis 중에 linking variable 의 index.
LinkingPerfidx = Info.LinkingPerfidx; % Linking Variable 이 몇 번째 index 인가?(FEA 에서)
LinkingInVar = Info.LinkingInVar; % Linking Variable 중 여기서 Variable
LinkingInPerf= Info.LinkingInPerf; % Linking Variable 중 여기서 Performance

Dspace = Info.Dspace;
lb = Dspace(1,:);
ub = Dspace(2,:);
Sample_in = Info.Sample_in; % Initial sampling using grid sampling in each subsystem.


%LinkingCDF = Info.LinkingCDF; % PATC 에서 사용할 Bottom-level 에서 전달될 CDF.
Seed = Info.Seed;
Obj_in = Info.Obj_in;
FEAFile = Info.FEAFile;
perfCritical = Info.perfCritical;
Objective = Info.Obj;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nc = 1.6; % size of local window
[S, Obj] = CBS(Sample_in,Obj_in,LinkingPerfidx,FEAFile,Var,nc,perfCritical);
ndv = size(S,2);
theta = ones(1,ndv); % Kriging parameter.

for i = 1 : size(Obj,2)
    [SurModel{i}, ~] = dacefit(S, Obj(:,i), @regpoly2, @corrgauss, theta);
end


if size(Obj_in,2) == 4
    d1 = 0:0.1:10;  d2 = 0:0.1:10;
    [D1,D2] = meshgrid(d1,d2);
    for n = 1:length(d1)
        for m = 1:length(d2)
            tmp = [D1(m,n),D2(m,n)];
            Z1_s(m,n) = predictor(tmp,SurModel{1}); % Draw Constraint 
            Z2_s(m,n) = predictor(tmp,SurModel{2}); % Draw Constraint 
            Z3_s(m,n) = predictor(tmp,SurModel{3}); % Draw Constraint 
        end
    end






    NumofStart = 25;
    tmpData = load('Det.mat');
    NumofDet = size(tmpData.SubsystemInfo{1,1}.Sample_in,1);
    ft=[0,0];
    figure(1)
    tmp_s = S;
    plot(tmp_s(1:NumofStart,1),tmp_s(1:NumofStart,2),'.','col','k','Markersize',20);
    hold on;
    plot(tmp_s(NumofStart+1:NumofDet,1),tmp_s(NumofStart+1:NumofDet,2),'x','col','r','Markersize',10,'Linewidth',1.5);
    hold on;
    try
        plot(tmp_s(NumofDet+1:end,1),tmp_s(NumofDet+1:end,2),'^','col','b','Markersize',10,'Linewidth',1.5);
        hold on;
    catch
    end
    plot(Var(1),Var(2),'s','col','k','Markersize',11,'linewidth',2);    
    contour(D1,D2,Z1_s,ft,':b','linewidth',1.5);
    hold on
    contour(D1,D2,Z3_s,ft,':r','linewidth',1.5);
    hold on
    contour(D1,D2,Z2_s,ft,':g','linewidth',1.5);
    hold on

  
    %legend('Const.1','Const.2','Const.3')
    xlabel('x_1');ylabel('x_2');
    Z1 = 1 -(D1.^2).*D2/20; % Draw Constraint 
    Z3 = 1 -80./(D1.^2 +8*D2 +5); % Draw Constraint
    Z4 = -1 + (0.9063*D1+0.4226*D2-6).^2+(0.9063*D1+0.4226*D2-6).^3-0.6*(0.9063*D1+0.4226*D2-6).^4-(-0.4226*D1+0.9063*D2); % Draw Constraint 2 
    Z5 = 1-(D1+D2-10).^2/30 -(D1-D2+10).^2/120; % f1(cost function)
    ft=[0 0];
    figure(1)
    contour(D1,D2,Z1,ft,'-b','linewidth',2);
    hold on
    contour(D1,D2,Z4,ft,'-r','linewidth',2);
    contour(D1,D2,Z3,ft,'-g','linewidth',2);
    legend('Initial sampling','Samples during ATC','Samples during PATC','Design point')
    hold off
    drawnow;
    print(['Figure',num2str(Globaliter)],'-dpng','-r600');
    
end


stepsize = 0.005;
Exitflag = -2;
while Exitflag==-2 || Exitflag==-1 || Exitflag==0
    stepsize = stepsize*0.7;
    options=optimoptions('fmincon','GradObj','off','GradConstr','on','TolCon',1e-3,'TolX',1e-3,'display','iter','algorithm','sqp','FiniteDifferenceStepSize',stepsize);
    [Var_opt,~,Exitflag,~]=fmincon(@ObjFile,Var,[],[],[],[],lb,ub,@ConstFile,options,Info,Objective,SurModel,v,w);
    Var = Var+Var.*(0.5*(rand(size(Var,1),size(Var,2))-0.5));
end

Var = Var_opt;




if ~isempty(LinkingPerfidx)
    
    for n = Localidx % Local variable 은 일단 Sampling.
        SmpVar(:,n) = norminv(Seed{n},Var(n),LocalVarStd(n)); % Local variable 들의 distribution 은 알고 있다고 가정. 일단은 모두 Normal(추후 변경 요망)
    end
    
    for n = Sharedidx % Local variable 은 일단 Sampling.
        SmpVar(:,n) = norminv(Seed{n},Var(n),SharedVarStd(n)); % Local variable 들의 distribution 은 알고 있다고 가정. 일단은 모두 Normal(추후 변경 요망)
    end
    
    if ~isempty(LinkingInVar) % Bottom level 로 부터 올라오는 임의의 CDF(kernel) 이 존재할 때.
        for n = 1 : length(LinkingInVar) % Bottom-level 에서 올라오는 임의의 CDF
            SmpCDF(:,n) = LinkingInVar(:,n)-mean(LinkingInVar(:,n))+Var(Coupledidx(n)); % 현재 Variable 값만큼 Shifting.
        end
        SmpVar = [SmpVar SmpCDF];
    end
         
    for i = 1 : length(LinkingPerfidx)
        LinkingOutNew(:,i) = predictor(SmpVar,SurModel{LinkingPerfidx(i)});
    end
    
else
    
    LinkingOutNew = [];
    
end

if ~isempty(Coupledidx)
    LinkingInNew = Var(Coupledidx);
else
    LinkingInNew = [];
end
if ~isempty(Sharedidx)   
    SharedVarNew = Var(Sharedidx);
else
    SharedVarNew = [];
end

VarNew = Var;
Sample_out= S;
Obj_out = Obj;

end

function Obj = ObjFile(Var,Info,Objective,SurModel,v,w)
    
    
    Obj = feval(Objective,Var); % Objective function 은 Mathematical 한 식으로 알 수 있다고 가정.
    
    RespPerf = [];
   
    
    Coupledidx = Info.Coupledidx; % Var 중에 어떤 것이 Coupled variable 인지
    Sharedidx = Info.Sharedidx; % Var 중에 어떤 것이 Shared variable 인지
    CoupledVar = Var(Coupledidx); % 여기서 Input 으로 작용하는 Coupling Variable
    SharedVar = Var(Sharedidx); % 여기서 Input 으로 작용하는 Shared Variable

    LinkingPerfidx = Info.LinkingPerfidx; % Linking Variable 이 몇 번째 index 인가?(FEA 에서)
    LinkingInVar = Info.LinkingInVar; % Linking Variable 중 여기서 Variable
    LinkingInPerf= Info.LinkingInPerf; % Linking Variable 중 여기서 Performance
    SharedIn = Info.SharedIn; % Linking Variable 중 여기서 Shared Variable
        
    CoupledLagidx = Info.CoupledLagidx;
    SharedLagidx = Info.SharedLagidx;
    LinkingLagidx = Info.LinkingLagidx;
    
        
    if ~isempty(LinkingPerfidx) % 여기서 Output 인 Coupling Variable
        
        for i = LinkingPerfidx
            RespPerf(i) = predictor(Var,SurModel{i});
        end
        
        RespPerf = RespPerf(LinkingPerfidx);
        
        Obj = Obj + v(LinkingLagidx)*(RespPerf-LinkingInPerf)' + (w(LinkingLagidx)).^2*(RespPerf-LinkingInPerf)'.^2;
        
    end
        
    
    if ~isempty(CoupledVar) % 여기서 Input 인 Coupling Variable
        
        for i = 1 : size(LinkingInVar,2)
            LinkingInVarMean(i) = mean(LinkingInVar(:,i)); % RBDO 이므로, Mean 으로 target 과 response 의 차이를 줄인다(Shifting)
        end
        Obj = Obj + v(CoupledLagidx)*(LinkingInVarMean-CoupledVar)' + (w(CoupledLagidx)).^2*(LinkingInVarMean-CoupledVar)'.^2;

        
        
    end
    
    if ~isempty(SharedVar) % 여기서 Input 인 Shared Variable
                
        Obj = Obj + v(SharedLagidx)*(SharedVar-SharedIn)' + (w(SharedLagidx)).^2*(SharedVar-SharedIn)'.^2;
        
        
    end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Reliabiliy Analysis  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [g,h,gradg,gradh] = ConstFile(Var,Info,Objective,SurModel,v,w)
    

    Localidx = Info.Localidx; % Var 중에 어떤 것이 Local variable 인지
    Coupledidx = Info.Coupledidx; % Var 중에 어떤 것이 Coupled variable 인지
    Sharedidx = Info.Sharedidx; % Var 중에 어떤 것이 Shared variable 인지
    CoupledVar = Var(Coupledidx); % 여기서 Input 으로 작용하는 Coupling Variable
    SharedVar = Var(Sharedidx); % 여기서 Input 으로 작용하는 Shared Variable

    LinkingPerfidx = Info.LinkingPerfidx; % Linking Variable 이 몇 번째 index 인가?(FEA 에서)
    LinkingInVar = Info.LinkingInVar; % Linking Variable 중 여기서 Variable
    LinkingInPerf= Info.LinkingInPerf; % Linking Variable 중 여기서 Performance
    SharedIn = Info.SharedIn; % Linking Variable 중 여기서 Shared Variable
    
    Seed = Info.Seed;
    BetaTarget = Info.Beta;
    LocalStd = Info.LocalVarStd;
    SharedStd = Info.SharedVarStd;
    perfCritical = Info.perfCritical;
    
    
    for n = Localidx % Local variable 은 일단 Sampling.
        SmpVar(:,n) = norminv(Seed{n},Var(n),LocalStd(n)); % Local variable 들의 distribution 은 알고 있다고 가정. 일단은 모두 Normal(추후 변경 요망)
    end
    
    for n = Sharedidx % Local variable 은 일단 Sampling.
        SmpVar(:,n) = norminv(Seed{n},Var(n),SharedStd(n)); % Local variable 들의 distribution 은 알고 있다고 가정. 일단은 모두 Normal(추후 변경 요망)
    end
    
    
    if ~isempty(LinkingInVar) % Bottom level 로 부터 올라오는 임의의 CDF(kernel) 이 존재할 때.
        for n = 1 : size(LinkingInVar,2) % Bottom-level 에서 올라오는 임의의 CDF
            SmpCDF(:,n) = LinkingInVar(:,n)-mean(LinkingInVar(:,n))+CoupledVar(n); % 현재 Variable 값만큼 Shifting.
        end
        SmpVar = [SmpVar SmpCDF];
    end
           
    sn = size(SmpVar,1);
    
    for i = 1 : length(SurModel) - length(LinkingPerfidx)
        Perf(:,i) = predictor(SmpVar,SurModel{i}); % 크리깅을 이용한 계산
        Pf = sum((Perf(:,i)-perfCritical(i))>0)/sn;
        g(i) = Pf - normcdf(-BetaTarget);
        
        for j = Localidx % Local variable 의 Sensitivity 계산
            gradg(j,i) = sum(((Perf(:,i)-perfCritical(i))>0).*((SmpVar(:,j)-Var(j))/LocalStd(j)^2))./sn;
        end
        
        for j = Sharedidx % Shared variable 의 Sensitivity 계산
            gradg(j,i) = sum(((Perf(:,i)-perfCritical(i))>0).*((SmpVar(:,j)-Var(j))/SharedStd(j)^2))./sn;
         end       
        
        for j = 1 : size(LinkingInVar,2) % Coupling Variable 의 Sensitivity 계산
            
            
            x_Data = SmpCDF(:,j);
            h=1.06*sqrt(var(x_Data))*sn^(-1/5); % Bandwidth
            sn_grad = 1000; % Sensitivity 계산하는 데 사용된 샘플 수 : 많지 않아도 정확하다. 너무 많으면 계산속도가 너무 느림.
            sn = size(x_Data,1); % 커널 만드는 데 사용된 샘플 수

            x1_kernelsmp = normrnd(0,h,[sn_grad,1]); % x1 sampling
            ker_idx = randi(sn,[sn_grad,1]); % shifting 하기 위한 임의의 정수 생성
            x1_kernelsmp = x1_kernelsmp + x_Data(ker_idx,:); % shifting MCsmp to each data
            x_cur(:,1) = x1_kernelsmp; % x1 from the kernel
          
           %% Analytical sensitivity
            x_kernel = @(x) sum(1/sn/h/sqrt(2*pi).*exp(-0.5*((x_Data-x)/h).^2)); % Unit kernel
            score_mean = zeros(sn_grad,1);
            for k = 1 : sn_grad
                result = x_kernel(x_cur(k,1)*ones(sn,1)); 
                score_mean(k) = sum(1./(h^3*result*sqrt(2*pi))/sn.*(x_cur(k,1).*ones(sn,1)-x_Data).*exp(-0.5*((x_cur(k,1).*ones(sn,1)-x_Data)/h).^2));
            end
            
            TmpSmpVar = SmpVar(ker_idx,:);
            TmpSmpVar(:,length(Localidx)+length(Sharedidx)+j) = x_cur(:,1);
            Perf_Tmp(:,i) = predictor(TmpSmpVar,SurModel{i}); % 크리깅을 이용한 계산  
            gradg(j+length(Localidx)+length(Sharedidx),i) = sum(((Perf_Tmp(:,i)-perfCritical(i))>0).*(score_mean))/sn_grad;

        end

    end

    h = [];
    gradh = [];

    
end
