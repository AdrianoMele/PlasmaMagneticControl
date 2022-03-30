function index = signalIndexByName(name,nameVector,indexVector)
% If name is a cell array of strings returns a vector of indexes
% If a name is present more than once in the nameVector only the first index will be returned
index = [];
if iscell(name)
  for i = 1:length(name)
    temp = indexVector(find(strcmp(nameVector,name{i})==1));
    if ~isempty(temp)
      index(end+1) = temp(1);
    end
  end
else
  temp = indexVector(find(strcmp(nameVector,name)==1));
  index = temp(1);
end
end