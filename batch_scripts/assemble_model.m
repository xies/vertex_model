function tisArr = assemble_model(DIR,max_time)
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

if nargin > 1,
    T = T(T < max_time);
end

n = numel(T);
% Preallocate
tisArr(1:n) = Tissue;

ind = 1;
for i = 1:10:n
    display(['Assembling time = ' num2str(T(i))]);
    tis.move_vts(Y(i,:), T(i));
    tisArr(ind) = Tissue(tis); % Make shallow copy of handle
    ind = ind + 1;
end

tisArr(isempty(tisArr)) = [];

save('-v7.3',[DIR '/model.mat'],'tisArr');
display(['Model assembled at: ' DIR '/model.mat']);

end