function Performance = FEA_subsystem2( Var )
    


    Const{1} =  @(x) 1 -(x(:,1).^2).*x(:,2)/20; % Constraint 1
    Const{2} =  @(x) 1 -(x(:,1)+x(:,2)-10).^2/30 -(x(:,1)-x(:,2)+10).^2/120; % Constraint 2
    
    Performance = [Const{1}(Var) Const{2}(Var)];

    
end

