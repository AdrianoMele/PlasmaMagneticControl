function idx = get_y_idx(y_type,y_name,flag_exact)
%idx = get_y_idx(y_type,y_name,[flag_exact])
%   Get index of the output y_name in the y_type cell array.

if nargin<3, flag_exact = 0; end
idx = contains(y_type(:,1),y_name);
idx = [y_type{idx,3}]'; % Group the indexes as a vector if more than one output has the same name

if flag_exact
  if ischar(y_name), y_name = {y_name}; end
  idx2 = zeros(numel(y_name),1);
  for i = 1:numel(y_name)
    ii = strcmp(y_type(idx,1),y_name{i});
    idx2(i) = idx(ii);    
  end
  idx = idx2;
end


