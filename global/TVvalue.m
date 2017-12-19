function value = TVvalue(patch)

value = sum(sum(abs(patch(:,1:end-1) - patch(:,2:end))));

end