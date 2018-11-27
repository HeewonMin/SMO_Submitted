% Subsystem 2
 
clear Info
Info.Num = 2;


Info.Localidx = [1]; % Design Variable �߿��� Local variable �� �� ��° index �� �ش��ϴ°� 

Info.Coupledidx = [2]; % Design Variable �߿��� Coupled variable �� �� ��° index �� �ش��ϴ°� 
Info.CoupledLagidx = [1]; % Lagrange multiplier �߿��� �� ��° index �� �ش��ϴ°�

Info.Sharedidx = []; % Design Variable �߿��� Shared variable �� �� ��° index �� �ش��ϴ°�
Info.SharedVarStd = []; % PATC ���� ���� Shared variable �� std ��.
Info.SharedVarType = {}; % PATC ���� ���� Shared variable �� ����Ÿ��.
Info.SharedLagidx = []; % Lagrange multiplier �߿��� �� ��° index �� �ش��ϴ°�

Info.LocalVar = [5]; % Design Variable �߿� Local variable.
Info.LocalVarStd = [0.5]; % PATC ���� ���� Local variable �� std ��.
Info.LocalVarType = {'Normal'}; % PATC ���� ���� Local variable �� ����Ÿ��.

Info.Var = [5 5]; % Variable �� initial value ��.

%% Output �κп��� Coupled Variable �κ�
Info.LinkingPerfidx = []; % FEA analysis �߿� linking variable �� index.
Info.LinkingLagidx = []; % Lagrange multiplier �߿��� �� ���� index �� �ش��ϴ°�?

Info.LinkingInVar = [5]; % �ٸ� subsystem ���κ��� ������ Linking Variable �� ���⼭�� Variable �� ���Ǵ� ��.
Info.LinkingInPerf = []; % �ٸ� subsystem ���κ��� ������ Linking Variable �� ���⼭�� Performance �� ���Ǵ� ��.
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
Info.Obj = @(x) x(1)+x(2); % �����Լ�; % �����Լ�.
Info.FEAFile = 'FEA_subsystem2'; % FEA �Լ�
Info.Obj_in = feval(Info.FEAFile,Info.Sample_in); % FEA �� �̿��ؼ� �������� �� Linking variable ���.
Info.perfCritical = zeros(1,size(Info.Obj_in,2)-length(Info.LinkingPerfidx)); % ���������� critical ��
SubsystemInfo{2} = Info; % ���⸦ subsystem �� ���ؼ� ���� �ٲ��ش�.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%