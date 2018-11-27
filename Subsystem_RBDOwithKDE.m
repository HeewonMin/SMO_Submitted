%%% General Subsystem in ATC %%%%
% Yongsu Jung, IDOLAB, Mechanical Engineering, KAIST
% E-mail : yongsu50@kaist.ac.kr
%% Input ����
% Var : ���� iteration �� design variable.
% Linking : ���� iteration �� Linking variable.
% ObjFile : �����Լ� function.
% ConstFile : �������� function.
% FEAFile : �������� �� Linking variable ����Լ�.
% v : Lagrange multiplier - ��� Lagrange multiplier �� �����ϹǷ� idx �� ������ index.
% w : Weight.
% Sample_in : ������ ���ø� ����.
% Obj_in : ������ ���ø� ����.(��������)
% Linking_in ; ������ ���ø� ����.(Coupling Variable)
% Bound : Design variable �� Upper and lower bound - Two row matrix.
% Options : fmincon �� ���� option.
% idx : ���� subsystem �� ���Ǵ� Lagrange multiplier �� ���� index.
%% Output ����
% LinkingNew : �̹� iteration �� ���� �� linkingvariable ���� - Distribution ���� �����ؾ� �Ѵ�.
% VarNew : �̹� iteration �� ���� �� design variable ����.
% Sample_out ; �̹� iteration �� ���� �� Sampling ����.
% LinkingIndex : Linking variable ���� Performance function �� �� ��°���� ��ġ�ϴ°�?
% SharedIndex : Linking variable �߿��� variable �� � variable ���� linking
% variable �ΰ�?


function [LinkingOutNew,LinkingInNew,SharedVarNew,VarNew,Sample_out,Obj_out] = Subsystem_RBDOwithKDE(Info,v,w)

Globaliter = Info.iter;

Localidx = Info.Localidx; % Var �߿� � ���� Local variable ����
Coupledidx = Info.Coupledidx; % Var �߿� � ���� Coupled variable ����
Sharedidx = Info.Sharedidx; % Var �߿� � ���� Shared variable ����

SharedVarStd = Info.SharedVarStd; % PATC ���� ���� Shared variable �� std ��.
SharedVarType = Info.SharedVarType; % PATC ���� ���� Shared variable �� ����

LocalVarStd = Info.LocalVarStd; % PATC ���� ���� local variable �� std ��.
LocalVarType = Info.LocalVarType; % PATC ���� ���� local variable �� ����

Var = Info.Var;

LinkingPerfidx = Info.LinkingPerfidx; % FEA analysis �߿� linking variable �� index.
LinkingPerfidx = Info.LinkingPerfidx; % Linking Variable �� �� ��° index �ΰ�?(FEA ����)
LinkingInVar = Info.LinkingInVar; % Linking Variable �� ���⼭ Variable
LinkingInPerf= Info.LinkingInPerf; % Linking Variable �� ���⼭ Performance

Dspace = Info.Dspace;
lb = Dspace(1,:);
ub = Dspace(2,:);
Sample_in = Info.Sample_in; % Initial sampling using grid sampling in each subsystem.


%LinkingCDF = Info.LinkingCDF; % PATC ���� ����� Bottom-level ���� ���޵� CDF.
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
    
    for n = Localidx % Local variable �� �ϴ� Sampling.
        SmpVar(:,n) = norminv(Seed{n},Var(n),LocalVarStd(n)); % Local variable ���� distribution �� �˰� �ִٰ� ����. �ϴ��� ��� Normal(���� ���� ���)
    end
    
    for n = Sharedidx % Local variable �� �ϴ� Sampling.
        SmpVar(:,n) = norminv(Seed{n},Var(n),SharedVarStd(n)); % Local variable ���� distribution �� �˰� �ִٰ� ����. �ϴ��� ��� Normal(���� ���� ���)
    end
    
    if ~isempty(LinkingInVar) % Bottom level �� ���� �ö���� ������ CDF(kernel) �� ������ ��.
        for n = 1 : length(LinkingInVar) % Bottom-level ���� �ö���� ������ CDF
            SmpCDF(:,n) = LinkingInVar(:,n)-mean(LinkingInVar(:,n))+Var(Coupledidx(n)); % ���� Variable ����ŭ Shifting.
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
    
    
    Obj = feval(Objective,Var); % Objective function �� Mathematical �� ������ �� �� �ִٰ� ����.
    
    RespPerf = [];
   
    
    Coupledidx = Info.Coupledidx; % Var �߿� � ���� Coupled variable ����
    Sharedidx = Info.Sharedidx; % Var �߿� � ���� Shared variable ����
    CoupledVar = Var(Coupledidx); % ���⼭ Input ���� �ۿ��ϴ� Coupling Variable
    SharedVar = Var(Sharedidx); % ���⼭ Input ���� �ۿ��ϴ� Shared Variable

    LinkingPerfidx = Info.LinkingPerfidx; % Linking Variable �� �� ��° index �ΰ�?(FEA ����)
    LinkingInVar = Info.LinkingInVar; % Linking Variable �� ���⼭ Variable
    LinkingInPerf= Info.LinkingInPerf; % Linking Variable �� ���⼭ Performance
    SharedIn = Info.SharedIn; % Linking Variable �� ���⼭ Shared Variable
        
    CoupledLagidx = Info.CoupledLagidx;
    SharedLagidx = Info.SharedLagidx;
    LinkingLagidx = Info.LinkingLagidx;
    
        
    if ~isempty(LinkingPerfidx) % ���⼭ Output �� Coupling Variable
        
        for i = LinkingPerfidx
            RespPerf(i) = predictor(Var,SurModel{i});
        end
        
        RespPerf = RespPerf(LinkingPerfidx);
        
        Obj = Obj + v(LinkingLagidx)*(RespPerf-LinkingInPerf)' + (w(LinkingLagidx)).^2*(RespPerf-LinkingInPerf)'.^2;
        
    end
        
    
    if ~isempty(CoupledVar) % ���⼭ Input �� Coupling Variable
        
        for i = 1 : size(LinkingInVar,2)
            LinkingInVarMean(i) = mean(LinkingInVar(:,i)); % RBDO �̹Ƿ�, Mean ���� target �� response �� ���̸� ���δ�(Shifting)
        end
        Obj = Obj + v(CoupledLagidx)*(LinkingInVarMean-CoupledVar)' + (w(CoupledLagidx)).^2*(LinkingInVarMean-CoupledVar)'.^2;

        
        
    end
    
    if ~isempty(SharedVar) % ���⼭ Input �� Shared Variable
                
        Obj = Obj + v(SharedLagidx)*(SharedVar-SharedIn)' + (w(SharedLagidx)).^2*(SharedVar-SharedIn)'.^2;
        
        
    end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Reliabiliy Analysis  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [g,h,gradg,gradh] = ConstFile(Var,Info,Objective,SurModel,v,w)
    

    Localidx = Info.Localidx; % Var �߿� � ���� Local variable ����
    Coupledidx = Info.Coupledidx; % Var �߿� � ���� Coupled variable ����
    Sharedidx = Info.Sharedidx; % Var �߿� � ���� Shared variable ����
    CoupledVar = Var(Coupledidx); % ���⼭ Input ���� �ۿ��ϴ� Coupling Variable
    SharedVar = Var(Sharedidx); % ���⼭ Input ���� �ۿ��ϴ� Shared Variable

    LinkingPerfidx = Info.LinkingPerfidx; % Linking Variable �� �� ��° index �ΰ�?(FEA ����)
    LinkingInVar = Info.LinkingInVar; % Linking Variable �� ���⼭ Variable
    LinkingInPerf= Info.LinkingInPerf; % Linking Variable �� ���⼭ Performance
    SharedIn = Info.SharedIn; % Linking Variable �� ���⼭ Shared Variable
    
    Seed = Info.Seed;
    BetaTarget = Info.Beta;
    LocalStd = Info.LocalVarStd;
    SharedStd = Info.SharedVarStd;
    perfCritical = Info.perfCritical;
    
    
    for n = Localidx % Local variable �� �ϴ� Sampling.
        SmpVar(:,n) = norminv(Seed{n},Var(n),LocalStd(n)); % Local variable ���� distribution �� �˰� �ִٰ� ����. �ϴ��� ��� Normal(���� ���� ���)
    end
    
    for n = Sharedidx % Local variable �� �ϴ� Sampling.
        SmpVar(:,n) = norminv(Seed{n},Var(n),SharedStd(n)); % Local variable ���� distribution �� �˰� �ִٰ� ����. �ϴ��� ��� Normal(���� ���� ���)
    end
    
    
    if ~isempty(LinkingInVar) % Bottom level �� ���� �ö���� ������ CDF(kernel) �� ������ ��.
        for n = 1 : size(LinkingInVar,2) % Bottom-level ���� �ö���� ������ CDF
            SmpCDF(:,n) = LinkingInVar(:,n)-mean(LinkingInVar(:,n))+CoupledVar(n); % ���� Variable ����ŭ Shifting.
        end
        SmpVar = [SmpVar SmpCDF];
    end
           
    sn = size(SmpVar,1);
    
    for i = 1 : length(SurModel) - length(LinkingPerfidx)
        Perf(:,i) = predictor(SmpVar,SurModel{i}); % ũ������ �̿��� ���
        Pf = sum((Perf(:,i)-perfCritical(i))>0)/sn;
        g(i) = Pf - normcdf(-BetaTarget);
        
        for j = Localidx % Local variable �� Sensitivity ���
            gradg(j,i) = sum(((Perf(:,i)-perfCritical(i))>0).*((SmpVar(:,j)-Var(j))/LocalStd(j)^2))./sn;
        end
        
        for j = Sharedidx % Shared variable �� Sensitivity ���
            gradg(j,i) = sum(((Perf(:,i)-perfCritical(i))>0).*((SmpVar(:,j)-Var(j))/SharedStd(j)^2))./sn;
         end       
        
        for j = 1 : size(LinkingInVar,2) % Coupling Variable �� Sensitivity ���
            
            
            x_Data = SmpCDF(:,j);
            h=1.06*sqrt(var(x_Data))*sn^(-1/5); % Bandwidth
            sn_grad = 1000; % Sensitivity ����ϴ� �� ���� ���� �� : ���� �ʾƵ� ��Ȯ�ϴ�. �ʹ� ������ ���ӵ��� �ʹ� ����.
            sn = size(x_Data,1); % Ŀ�� ����� �� ���� ���� ��

            x1_kernelsmp = normrnd(0,h,[sn_grad,1]); % x1 sampling
            ker_idx = randi(sn,[sn_grad,1]); % shifting �ϱ� ���� ������ ���� ����
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
            Perf_Tmp(:,i) = predictor(TmpSmpVar,SurModel{i}); % ũ������ �̿��� ���  
            gradg(j+length(Localidx)+length(Sharedidx),i) = sum(((Perf_Tmp(:,i)-perfCritical(i))>0).*(score_mean))/sn_grad;

        end

    end

    h = [];
    gradh = [];

    
end
