function Performance = FEA_subsystem1( Var )
    

    Const{1} = @(x) 1-x(:,1).^2.*x(:,2)./20; % Constraint 1
    Const{2} = @(x) 1-(80./(x(:,1).^2+8*x(:,2)+5)); % Constraint 2
    Const{3} = @(x) -1 + (0.9063*x(:,1)+0.4226*x(:,2)-6).^2+(0.9063*x(:,1)+0.4226*x(:,2)-6).^3-0.6*(0.9063*x(:,1)+0.4226*x(:,2)-6).^4-(-0.4226*x(:,1)+0.9063*x(:,2));
    LinkingConst{1} = @(x) (2*x(:,1).^2-x(:,2))./10; % Linking variable
    Performance = [Const{1}(Var) Const{2}(Var) Const{3}(Var)];
    Performance = [Performance LinkingConst{1}(Var)];
    
end

