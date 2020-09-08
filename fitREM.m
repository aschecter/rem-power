function Output = fitREM(data,x0)

tic

%-------------------------------------------------------------------------%
% Initialize
%-------------------------------------------------------------------------%

numStats = length(x0);
numEvents = size(data,1);

statistics = sparse(data(:,2:end-1));
events = data(:,1);
event_idx = data(:,end);
clear data

theta = x0;

%-------------------------------------------------------------------------%
% Fit Sequence
%-------------------------------------------------------------------------%

%%% Run MATLAB's prebuilt solver for unconstrained non-linear problems
% options = optimset('GradObj','on','Hessian','on','Display','off','TolX',10^-6,'TolFun',10^-6,'MaxIter',500);
options = optimoptions('fminunc','Algorithm','trust-region','Display','iter','SpecifyObjectiveGradient',true,'HessianFcn','objective');
[x,~,exitflag,output,~,~] = fminunc(@LGH,theta,options);

elapsedTime = toc;
%-------------------------------------------------------------------------%
% Save Output
%-------------------------------------------------------------------------%
N = length(x);
Output.parameters = x;     

[fval,grad,hess] = LGH(x);
Output.grad = grad;
Output.hess = hess;
% save the final likelihood value (negative since the algorithm solved a
% minimization problem)    
Output.likelihood = -fval;

% Optimization information       
Output.status = exitflag;
Output.info = output;

% SD and Var
var = inv(Output.hess);
std = zeros(N,1);
for ss = 1:N
    if var(ss,ss)>=0
        std(ss) = sqrt(var(ss,ss));
    end
end
Output.stds = std;
Output.var = var;

% p-values
Output.pval = 2*(1-normcdf(abs(Output.parameters./Output.stds)));

% Compute Correlations
corr = eye(N);
for i = 1:(N-1)
    for j = (i+1):N
        if (std(i)*std(j) ~= 0) && isfinite(var(i,j))
            corr(i,j) = var(i,j)/(std(i)*std(j));
            corr(j,i) = corr(i,j);
        end
    end
end
Output.correlation = corr;

% AIC and BIC
AIC = 2*N - 2*Output.likelihood;
BIC = -2*Output.likelihood + N*log(numEvents);
Output.AIC = AIC;
Output.BIC = BIC;

% Time performance
Output.time = elapsedTime;

function [L,gradient,hessian] = LGH(theta)
    gradient = zeros(numStats,1);
    hessian = zeros(numStats);
    [~,~,c] = unique(event_idx);
    
    % Rate Matrix
    R = zeros(size(statistics,1),1);
    for k = 1:numStats
        R = R + theta(k)*statistics(:,k);    % combine parameter and statistic values
    end
    rate = exp(R);                       % Take the exponential

    % Update loglikelihood function 
    L = sum(R(events==1)) - sum(log(accumarray(c,rate)));
    clear R

    % compute gradient                
    for k = 1:numStats         
        gradient(k) = sum(statistics(events==1,k)) - sum(accumarray(c,statistics(:,k).*rate)./accumarray(c,rate));
    end 

    % compute hessian
    for xh = 1:numStats
        for y = xh:numStats
            hessian(xh,y) = -sum(accumarray(c,(statistics(:,xh).*statistics(:,y)).*rate)./accumarray(c,rate)) + sum((accumarray(c,statistics(:,xh).*rate).*accumarray(c,statistics(:,y).*rate))./(accumarray(c,rate).^2));
            hessian(y,xh) = hessian(xh,y);
        end
    end 
    clear rate

    % negate the values since it is a maximization problem
    L = -L;
    gradient = -gradient;
    hessian = -hessian;
        

end

end