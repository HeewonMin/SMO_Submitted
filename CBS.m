function [Samples, Obj] = CBS(S,Obj,LinkingPerf,FEAFile,x,nc,perfCritical)


NumofObj = size(Obj,2); % The number of performance.
ndv = size(S,2); % The number of variables.
numofinitial = length(S(:,1)); % current number of samples
minbound = x-nc; % size of local window
maxbound = x+nc; % size of local window
numpredic = round((4*10^4)^(1/ndv)); % Grid sampling 
theta = ones(1,ndv); % Kriging parameter.
lob = 0.1.*ones(1,ndv); % Kriging parameter.
upb = 20.*ones(1,ndv); % Kriging parameter.
Iter_limit = 100; % CBS maximum iteration
X_Gridsamp = gridsamp([minbound;maxbound],numpredic); 
tmp_d = X_Gridsamp-repmat(x,size(X_Gridsamp,1),1);
dist = sqrt(sum(tmp_d.^2,2));
k = find(dist<nc);
X_Gridsamp = X_Gridsamp(k,:); 

for iter = 1:Iter_limit
    
    for i = 1 : NumofObj-length(LinkingPerf) 
        [dmodel{i},~]=dacefit(S,Obj(:,i),@regpoly2,@corrgauss,theta,lob,upb); 
        [Obj_X(:,i),mse(:,i)] = predictor(X_Gridsamp, dmodel{i});
        mean_mse(i) = abs(max(mse(:,i))/mean(Obj(1:end,i))); % accuracy indicator
        
 
    end
    accuracy = max(mean_mse);
    if iter == 1 && accuracy < 10^(-6) 
        Samples = S;
        return;
    end
    
    for i = 1:size(X_Gridsamp,1)
        dd = S-repmat(X_Gridsamp(i,:),size(S,1),1);
        summation = zeros(size(S,1),1);
        for j= 1:ndv
            summation = summation + dd(:,j).^2;
        end
        distance(i) = min(summation.^0.5 ); 

        tmp = 0; 
        
        if Obj_X(i,:) <= perfCritical
            for k = 1:NumofObj-length(LinkingPerf)
                tmp = tmp + normpdf((perfCritical(k)-(Obj_X(i,k)))/sqrt(mse(i,k)));
            end      
        else
             tmp = 0;
        end
        
        CBS_criterion(i) = tmp; 
        CBS(i) = distance(i)*CBS_criterion(i); 

    end    

    [~, index1] = max(CBS);
    if all(CBS==0)
        Samples = S;
        return;
    end
    newpoint = X_Gridsamp(index1,:); 

    if any(S == repmat(newpoint,size(S,1),1)) 
        CBS(index1) = [];
        [~, index2]=max(CBS); 
        newpoint2 = X_Gridsamp(index2,:); 
        if newpoint == newpoint2
            Samples = S;
            return;
        end
    end
       
    S = [S ; newpoint];
    Obj(numofinitial+iter,:) = feval(FEAFile,newpoint);
    
    for i = 1 : NumofObj-length(LinkingPerf)
        error_CBS(iter,i) = abs((predictor(S(end,:),dmodel{i})-Obj(numofinitial+iter,i)))./(max(Obj_X(:,i))-min(Obj_X(:,i)));
    end

    
    if max(error_CBS(iter,:)) < 1e-2
        break;
    end
    
end

Samples = S;

end

