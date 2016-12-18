function ranked = rankavec(vector)
% return the rank of each element in a vector

ranked = zeros(length(vector),1);

for i=1:length(vector),
        ranked(i) = length(find(vector<vector(i)))+1;
end


