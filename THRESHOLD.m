%%%%%%%%%%%%%%%%%%%%%%%%%%%% THRESHOLD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ phiproj, etha] = THRESHOLD(phi, beta) % ͶӰ
%% VOLUME-PRESERVING HEAVISIDE PROJECTION
[length_y,length_x] = size(phi);
phi(phi<0) = 0;phi(phi>1) = 1;
phi = phi(:); eta = [0,1];
while eta(2) - eta(1) > 0.0001 
    phiproj = zeros(1,length(phi));
    etha = eta(1);
    for i = 1 : length(phi)
        if phi(i) < etha
            phiproj(i) = etha*(exp(-beta*(1-phi(i)/etha))....
                - (1-phi(i)/etha)*exp(-beta));
        else
            phiproj(i) = (1-etha)*(1-exp(-beta*(phi(i)-etha)/(1-etha))....
                + (phi(i)-etha)*exp(-beta)/(1-etha))+etha;
        end
    end
    error1 = sum(sum(phi)) - sum(sum(phiproj));
 
    phiproj = zeros(1,length(phi));
    etha = (eta(1) + eta(2))/2;
    for i = 1 : length(phi)
        if phi(i) < etha
            phiproj(i) = etha*(exp(-beta*(1-phi(i)/etha))....
                - (1-phi(i)/etha)*exp(-beta));
        else
            phiproj(i) = (1-etha)*(1-exp(-beta*(phi(i)-etha)/(1-etha))....
                + (phi(i)-etha)*exp(-beta)/(1-etha))+etha;
        end
    end
    error2 = sum(sum(phi)) - sum(sum(phiproj));
    
    if error1*error2 <= 0
        eta(2) = (eta(1)+eta(2))/2;
    else
        eta(1) = (eta(1)+eta(2))/2;
    end
end
phiproj = zeros(1,length(phi));
etha = (eta(1)+eta(2))/2;
for i = 1 : length(phi)
    if phi(i) < etha
        phiproj(i) = etha*(exp(-beta*(1-phi(i)/etha))....
            - (1-phi(i)/etha)*exp(-beta));
    else
        phiproj(i) = (1-etha)*(1-exp(-beta*(phi(i)-etha)/(1-etha))....
            + (phi(i)-etha)*exp(-beta)/(1-etha))+etha;
    end
end
phiproj = reshape(phiproj,length_y,length_x);
end