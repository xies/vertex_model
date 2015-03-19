function t = tension(forces)
% given n forces, find the "tension" using the components

Fx = forces(:, 1);
Fy = forces(:, 2);

Fx_pos = Fx(Fx > 0);
Fy_pos = Fy(Fy > 0);

Fx_final_pos = sum(Fx_pos);
Fy_final_pos = sum(Fy_pos);

t = sqrt(Fx_final_pos^2 + Fy_final_pos^2);

% Fx_neg = Fx(Fx < 0);
% Fy_neg = Fy(Fy < 0);
% 
% Fx_final_neg = sum(Fx_neg);
% Fy_final_neg = sum(Fy_neg);
% 
% if Fx_final_pos ~= Fx_final_neg
%     disp('bad');
% end
% 
% if Fy_final_pos ~= Fy_final_neg
%     disp('bad');
% end   
%     