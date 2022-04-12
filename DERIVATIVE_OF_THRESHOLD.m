%%%%%%%%%%%%%%%%%%%%%%%% DERIVATIVE_OF_THRESHOLD %%%%%%%%%%%%%%%%%%%%%%%%%%
function [ dproj ] = DERIVATIVE_OF_THRESHOLD( phi, beta, etha) % 投影导数
[length_y,length_x] = size(phi);
phi = phi(:);
dproj = zeros(1,length(phi));
for i = 1 : length(phi)
    if phi(i) < etha
        dproj(i) = beta * exp(-beta*(1-phi(i)/etha))....
            + exp(-beta);
    else
        dproj(i) = beta * exp(-beta*(phi(i)-etha)/(1-etha))....
            + exp(-beta);
    end
end
dproj = reshape(dproj,length_y,length_x);
dproj(phi<0) = 0;dproj(phi>1) = 0;
end