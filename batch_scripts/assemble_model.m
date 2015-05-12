function tisArr = assemble_model(DIR)
%ASSEMBLE_MODEL
% Loads all .mat files (of Tissues) within a single directory and assemble
% into a tissue array. And then sort tissueArray by time stamp.
%
% Usage: tissueArray = assemble_model(DIR)

% Grab directory contents
s = what(DIR);
matfiles = s.mat;
% Preallocate
n = numel(matfiles);
tisArr(1:n) = Tissue;

for i = 1:n
    load( [DIR '/' matfiles{i}] );
    tisArr(i) = tis;
end

tisArr = tisArr.sort;

end