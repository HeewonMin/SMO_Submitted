%%% General Subsystem in ATC %%%%
% Yongsu Jung, IDOLAB, Mechanical Engineering, KAIST
% E-mail : yongsu50@kaist.ac.kr
%% Input ����
% Coupledidx : Variable �߿� � ���� Coupled Variable ������ ���� �ε���
% Sharedidx : Variable �߿� � ���� Shared Variable ������ ���� �ε���
% LinkingPerfidx : FEA ���� ������ Performance �߿� �� ��° �ε����� Coupling Varible �ΰ�?

% LinkingInVar : Linking Variable �� ���⼭ Variable (�ٸ��ý��ۿ��� ���� ��)
% LinkingInPerf : Linking Variable �� ���⼭ Performance (�ٸ��ý��ۿ��� ���� ��)
% SharedIn  : Linking Variable �� ���⼭ Shared Variable (�ٸ��ý��ۿ��� ���� ��)

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
%% Output ����



function [LinkingOutNew,LinkingInNew,SharedVarNew,VarNew,Sample_out,Obj_out] = Subsystem(Info,v,w)


  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Setting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Globaliter = Info.iter;

Coupledidx = Info.Coupledidx; % Var �߿� � ���� Coupled variable ����
Sharedidx = Info.Sharedidx; % Var �߿� � ���� Shared variable ����

SharedVarStd = Info.SharedVarStd; % PATC ���� ���� Shared variable �� std ��.
SharedVarType = Info.SharedVarType; % PATC ���� ���� Shared variable �� ����

LocalVarStd = Info.LocalVarStd; % PATC ���� ���� local variable �� std ��.
LocalVarType = Info.LocalVarType; % PATC ���� ���� local variable �� ����

Var = Info.Var;

LinkingPerfidx = Info.LinkingPerfidx; % FEA analysis �߿� linking variable �� index.


Dspace = Info.Dspace;
lb = Dspace(1,:);
ub = Dspace(2,:);
Sample_in = Info.Sample_in; % Initial sampling using grid sampling in each subsystem.
options = Info.options;

Beta = Info.Beta;
Obj_in = Info.Obj_in;
FEAFile = Info.FEAFile;
perfCritical = Info.perfCritical;
Objective = Info.Obj;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nc = 0.3; % size of local window
S = Sample_in;
Obj = Obj_in;
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
    ft=[0,0];
    figure(1)
    tmp_s = S;
    plot(tmp_s(1:NumofStart,1),tmp_s(1:NumofStart,2),'.','col','k','Markersize',20);
    hold on;
    plot(tmp_s(NumofStart+1:end,1),tmp_s(NumofStart+1:end,2),'x','col','r','Markersize',10,'Linewidth',1.5);
    hold on;
    plot(Var(1),Var(2),'s','col','k','Markersize',11,'linewidth',2);
    hold on;
    contour(D1,D2,Z1_s,ft,':b','linewidth',1.5);
    hold on
    contour(D1,D2,Z3_s,ft,':r','linewidth',1.5);
    hold on
    contour(D1,D2,Z2_s,ft,':g','linewidth',1.5);
    hold on

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


    legend('Initial sampling','Samples during ATC','Design point')
    %legend('Const.1(Kriging)','Const.2(Kriging)','Const.3(Kriging)','Const.1(True)','Const.2(True)','Const.3(True)','Initial Sampling','Sampling during deterministic ATC')
    hold off
    drawnow;
    print(['Figure',num2str(Globaliter),'Deter'],'-dpng','-r600');
end







Var=fmincon(@ObjFile,Var,[],[],[],[],lb,ub,@ConstFile,options,Info,Objective,SurModel,v,w);
[S, Obj] = CBS(Sample_in,Obj_in,LinkingPerfidx,FEAFile,Var,nc,perfCritical);



if ~isempty(LinkingPerfidx)
    for i = 1 : length(LinkingPerfidx)
        LinkingOutNew(i) = predictor(Var,SurModel{LinkingPerfidx(i)});
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
        
        Obj = Obj + v(LinkingLagidx).*(RespPerf-LinkingInPerf) + (w(LinkingLagidx)).^2.*(RespPerf-LinkingInPerf).^2;
        
    end
        
    
    if ~isempty(CoupledVar) % ���⼭ Input �� Coupling Variable
        
        
        Obj = Obj + v(CoupledLagidx).*(LinkingInVar -CoupledVar) + (w(CoupledLagidx)).^2.*(LinkingInVar -CoupledVar).^2;

        
        
    end
    
    if ~isempty(SharedVar) % ���⼭ Input �� Shared Variable
                
        Obj = Obj + v(SharedLagidx).*(SharedVar-SharedIn) + (w(SharedLagidx)).^2.*(SharedVar-SharedIn).^2;
        
        
    end
    


end
function [g,h] = ConstFile(Var,Info,Objective,SurModel,v,w)

    perfCritical = Info.perfCritical;
    for i = 1 : length(SurModel)
        g(i) = predictor(Var,SurModel{i}); % ũ������ �̿��� ���
    end
    LinkingPerfidx = Info.LinkingPerfidx;
    g(LinkingPerfidx) = []; % �������Ǹ� ����� ���� Linking variable �� ����.
    g = g-perfCritical;
    h = [];
    
end
