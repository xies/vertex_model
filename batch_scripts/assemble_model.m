function tisArr = assemble_model(DIR)
%ASSEMBLE_MODEL
% Loads all .mat files (of Tissues) within a single directory and assemble
% into a tissue array. And then sort tissueArray by time stamp.
%
% Usage: tissueArray = assemble_model(DIR)

% Grab directory contents
T = csvread([DIR '/times.csv']);
Y = csvread([DIR '/vertices.csv']);
tis = load([DIR '/model_t_0.mat']);
tis = tis.tis;

n = numel(T);
% Preallocate
tisArr(1:n) = Tissue;

for i = 1:n
    tis.evolve(Y(i,:), T(i));
    tisArr(i) = Tissue(tis); % Make shallow copy
end

save([DIR 'model.mat'],'tisArr');

end