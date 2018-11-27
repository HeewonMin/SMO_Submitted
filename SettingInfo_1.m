% Subsystem 1
 
clear Info
Info.Num = 1;

%% Input 
Info.Localidx = [1 2]; % Design Variable 중에서 Local variable 이 몇 번째 index 에 해당하는가 

Info.Coupledidx = []; % Design Variable 중에서 Coupled variable 이 몇 번째 index 에 해당하는가 
Info.CoupledLagidx = []; % Lagrange multiplier 중에서 몇 번째 index 에 해당하는가

Info.Sharedidx = []; % Design Variable 중에서 Shared variable 이 몇 번째 index 에 해당하는가
Info.SharedVarStd = []; % PATC 에서 쓰일 Shared variable 의 std 값.
Info.SharedVarType = {}; % PATC 에서 쓰일 Shared variable 의 분포타입.
Info.SharedLagidx = []; % Lagrange multiplier 중에서 몇 번째 index 에 해당하는가

Info.LocalVar = [5 5]; % Design Variable 중에 Local variable.
Info.LocalVarStd = [0.5 0.5]; % PATC 에서 쓰일 Local variable 의 std 값.
Info.LocalVarType = {'Normal','Normal'}; % PATC 에서 쓰일 Local variable 의 분포타입.

Info.Var = [5 5]; % Variable 의 initial value 값.

%% Output 부분에서 Coupled Variable 부분
Info.LinkingPerfidx = [4]; % FEA analysis 중에 linking variable 의 index.
Info.LinkingLagidx = [1]; % Lagrange multiplier 중에서 몇 번재 index 에 해당하는가?

Info.LinkingInVar = []; % 다른 subsystem 으로부터 들어오는 Linking Variable 중 여기서는 Variable 로 사용되는 값.
Info.LinkingInPerf = [5]; % 다른 subsystem 으로부터 들어오는 Linking Variable 중 여기서는 Performance 로 사용되는 값.
Info.SharedIn = []; % Linking Variable 중 여기서 Shared Variable 로 사용되는 값.

%% Surrogate modeling 을 위한 Initial Sampling
Info.Dspace= [[0 0];[10 10]];
InitialSampling = 5;
Info.Sample_in=gridsamp(Info.Dspace,InitialSampling); % Initial sampling using grid sampling in each subsystem.


%% 최적화 돌리기 위한 Option
Info.options = options;


%% PATC 를 위한 CDF 구성 및 Seed
Pf = 0.05; % Probability of Failure.
Info.Beta=-norminv(Pf); % Target index.

%% Objective function 과 FEA 설정(제한조건 및 Linking variable 계산)
Info.Obj = @(x) -(x(:,1)+x(:,2)-10).^2/30-(x(:,1)-x(:,2)+10).^2./120; % 목적함수.
Info.FEAFile = 'FEA_subsystem1'; % FEA 함수
Info.Obj_in = feval(Info.FEAFile,Info.Sample_in); % FEA 를 이용해서 제한조건 및 Linking variable 계산.
Info.perfCritical = zeros(1,size(Info.Obj_in,2)-length(Info.LinkingPerfidx)); % 제한조건의 critical 값
SubsystemInfo{1} = Info; % 여기를 subsystem 에 대해서 숫자 바꿔준다.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%