function sequence = createSequence4(riskSet,events,effect,effectSize)

%
% riskSet is a number corresponding to number of individuals (assume square
% matrix). events is number of events in sequence. effect is an index from
% 1:5 indicating which statistic has a non-zero effect. effectSize is the
% strength of effect (given as exp(theta))
%


%%% Initialization
sequence = zeros(events,2);  % empty sequence vector
weight = zeros(riskSet);     % weight matrix

idx = randsample(length(weight(:)),1,'true');  % first event
[row,col] = ind2sub(size(weight),idx);
weight(row,col) = 1;
sequence(1,1) = row;
sequence(1,2) = col;

%%% Set parameters
theta = zeros(5,1);
theta(effect) = log(effectSize);

%%% Generate sequence

for e = 2:events

    % Create Stats
    stats1 = weight;
    stats2 = weight';
    stats3 = repmat(sum(weight,2),1,riskSet) - diag(diag(repmat(sum(weight,2),1,riskSet)));
    stats4 = repmat(sum(weight,1),riskSet,1) - diag(diag(repmat(sum(weight,1),riskSet,1)));
    stats5 = 2*sqrt(weight*weight' - diag(diag(weight*weight')));
    
  
    % Create rates
    rate = (theta(1)*stats1 + theta(2)*stats2 + theta(3)*stats3 + theta(4)*stats4 + theta(5)*stats5);
    rate = exp(rate);
    rate = rate - diag(diag(rate));  % no self loop
    
    % Draw sequence
    if min(min(rate))==0 && max(max(rate))==0
        rate = ones(riskSet) - eye(riskSet);    % if all rates are zero (error) make all rates one to make dyads equally likely
    end
    I = randsample(length(rate(:)),1,'true',rate(:));   % randomly select index
    [row,col] = ind2sub(size(rate),I);                  % convert to row/column indices
    
    % Set sender and receiver
    sequence(e,1) = row;
    sequence(e,2) = col;
    
    % Update weight matrix
    weight(row,col) = weight(row,col) + 1;    
    
end