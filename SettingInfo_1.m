% Subsystem 1
 
clear Info
Info.Num = 1;

%% Input 
Info.Localidx = [1 2]; % Design Variable �߿��� Local variable �� �� ��° index �� �ش��ϴ°� 

Info.Coupledidx = []; % Design Variable �߿��� Coupled variable �� �� ��° index �� �ش��ϴ°� 
Info.CoupledLagidx = []; % Lagrange multiplier �߿��� �� ��° index �� �ش��ϴ°�

Info.Sharedidx = []; % Design Variable �߿��� Shared variable �� �� ��° index �� �ش��ϴ°�
Info.SharedVarStd = []; % PATC ���� ���� Shared variable �� std ��.
Info.SharedVarType = {}; % PATC ���� ���� Shared variable �� ����Ÿ��.
Info.SharedLagidx = []; % Lagrange multiplier �߿��� �� ��° index �� �ش��ϴ°�

Info.LocalVar = [5 5]; % Design Variable �߿� Local variable.
Info.LocalVarStd = [0.5 0.5]; % PATC ���� ���� Local variable �� std ��.
Info.LocalVarType = {'Normal','Normal'}; % PATC ���� ���� Local variable �� ����Ÿ��.

Info.Var = [5 5]; % Variable �� initial value ��.

%% Output �κп��� Coupled Variable �κ�
Info.LinkingPerfidx = [4]; % FEA analysis �߿� linking variable �� index.
Info.LinkingLagidx = [1]; % Lagrange multiplier �߿��� �� ���� index �� �ش��ϴ°�?

Info.LinkingInVar = []; % �ٸ� subsystem ���κ��� ������ Linking Variable �� ���⼭�� Variable �� ���Ǵ� ��.
Info.LinkingInPerf = [5]; % �ٸ� subsystem ���κ��� ������ Linking Variable �� ���⼭�� Performance �� ���Ǵ� ��.
Info.SharedIn = []; % Linking Variable �� ���⼭ Shared Variable �� ���Ǵ� ��.

%% Surrogate modeling �� ���� Initial Sampling
Info.Dspace= [[0 0];[10 10]];
InitialSampling = 5;
Info.Sample_in=gridsamp(Info.Dspace,InitialSampling); % Initial sampling using grid sampling in each subsystem.


%% ����ȭ ������ ���� Option
Info.options = options;


%% PATC �� ���� CDF ���� �� Seed
Pf = 0.05; % Probability of Failure.
Info.Beta=-norminv(Pf); % Target index.

%% Objective function �� FEA ����(�������� �� Linking variable ���)
Info.Obj = @(x) -(x(:,1)+x(:,2)-10).^2/30-(x(:,1)-x(:,2)+10).^2./120; % �����Լ�.
Info.FEAFile = 'FEA_subsystem1'; % FEA �Լ�
Info.Obj_in = feval(Info.FEAFile,Info.Sample_in); % FEA �� �̿��ؼ� �������� �� Linking variable ���.
Info.perfCritical = zeros(1,size(Info.Obj_in,2)-length(Info.LinkingPerfidx)); % ���������� critical ��
SubsystemInfo{1} = Info; % ���⸦ subsystem �� ���ؼ� ���� �ٲ��ش�.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%