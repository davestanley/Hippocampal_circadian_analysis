function h = plotIM (traw,xraw,IM)
% Part of a series of functions designed to operate on index matrix, IM.
% IM is a matrix whereby each row is an index vector, used to index input
% data (traw, xraw)
%

% Make sure IM are aligned as rows
if isvector(IM)
    IM = IM(:)';
end

h = cell(1,size(IM,1));

colourarr = get_colourarr;

for i = 1:size(IM,1)
    ind = IM(i,:);
    hold on; h{i} = plot(traw(ind),xraw(ind),colourarr(i));
end


end