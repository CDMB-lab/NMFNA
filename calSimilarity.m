% Algorithm for calculating the similarity of community.
function [ Sim ] = calSimilarity( comunities )
communityCounts = length(comunities);
Sim = zeros(communityCounts);
for i = 1:(communityCounts-1)
    community_i = comunities{i};
    for j = (i+1):communityCounts
        community_j = comunities{j};
        if ((~isempty(community_i)) && (~isempty(community_j)))
            Sim(i, j) = length(intersect(community_i,community_j))/min(length(community_i),length(community_j));
        end
    end
end
end