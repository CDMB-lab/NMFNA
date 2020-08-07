% Algorithm for Assigning the community members and truing the community to obtain more accurate results
function [ communityFinal,mSim ] = communitySelection( G1, xita )
[N, K] = size(G1);
S_mean = mean(G1);
S_std = std(G1, 0, 1);

candidatecommmunity = cell(K, 1);
communitySignal = zeros(K, 1);

for k = 1:K
    candidatecommmunity{k} = find(G1(:, k) > S_mean(k) + xita*S_std(k)); % Z-score>=t
    communitySignal(k) = mean(G1(candidatecommmunity{k}, k));
end

% Calculating the similarity of community
Sim = calSimilarity( candidatecommmunity );
communityFinal = candidatecommmunity;

% Obtaining and merging the communities
for i = 1:size(Sim, 1)-1
    for j = (i+1):size(Sim, 2)
        if Sim(i,j)>0.5
            [Y, I] = max([communitySignal(i), communitySignal(j)]);
            if I == 1
                communityFinal{j} = [];
                communitySignal(j) = 0;
                Sim(j, :) = zeros(1, size(Sim, 2));
                Sim(:, j) = zeros(size(Sim, 1), 1);
            else
                communityFinal{i} = [];
                communitySignal(i) = 0;
                Sim(i, :) = zeros(1, size(Sim, 2));
                Sim(:, i) = zeros(size(Sim, 1), 1);
            end
        end
    end
end

% Keeping the communities with no less than 5 node
i = 1;
while i ~= length(communityFinal)+1
    if isempty(communityFinal{i}) || (length(communityFinal{i})<5)
        communityFinal(i) = [];
        i = i - 1;
    end
    i = i + 1;
end

mSim=sum((sum(Sim))')/sum(sum(Sim~=0));
end