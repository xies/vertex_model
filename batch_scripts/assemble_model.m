function tisArr = assemble_model(DIR)
% Loads all .mat files (of Tissues) within a single directory and assemble
% into a tissue array. And then sort tissueArray by time stamp.
%
% Usage: tissueArray = assemble_model(DIR)

% Grab directory contents
s = what(DIR);
tisArr(1:n) = Tissue;

for i = 1:n
    load([DIR '/model_step_' num2str(i) '.mat'])
    tisArr(i) = tis;
end

end