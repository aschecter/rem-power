function data = createData2(sequence,riskSet)

events = size(sequence,1);
weight = zeros(riskSet,riskSet);
data = [];
s = sequence(1,1);   % sender
r = sequence(1,2);   % receiver
idx = eye(riskSet);

% Update Weight
weight(s,r) = 1;

    
for e = 2:events
     
    W = weight;
    
    s = sequence(e,1);   % sender
    r = sequence(e,2);   % receiver 
    
    % Create Stats
    stats1 = W;
    stats2 = W';
    stats3 = repmat(sum(W,2),1,riskSet) - diag(diag(repmat(sum(W,2),1,riskSet)));
    stats4 = repmat(sum(W,1),riskSet,1) - diag(diag(repmat(sum(W,1),riskSet,1)));
    stats5 = 2*sqrt(W*W' - diag(diag(W*W')));
    
%     % Standardize
%     stats1 = stats1/(e-1);
%     stats2 = stats2/(e-1);
%     stats3 = stats3/(e-1);
%     stats4 = stats4/(e-1);
%     stats5 = stats5/(e-1);

    temp = zeros(riskSet);
    temp(s,r) = 1;
    data = [data;temp(~idx) stats1(~idx) stats2(~idx) stats3(~idx) stats4(~idx) stats5(~idx) (e-1)*ones(length(temp(~idx)),1)];    
    
    % Update Weight    
    weight(s,r) = weight(s,r) + 1;
       
end

