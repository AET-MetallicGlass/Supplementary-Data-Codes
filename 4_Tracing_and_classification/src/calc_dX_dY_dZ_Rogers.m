function [dX, dY, dZ] = calc_dX_dY_dZ_Rogers(Orders,Coeff)
    P_200 = Coeff(Orders(:,1)==2 & Orders(:,2)==0 & Orders(:,3)==0);
    P_020 = Coeff(Orders(:,1)==0 & Orders(:,2)==2 & Orders(:,3)==0);
    P_002 = Coeff(Orders(:,1)==0 & Orders(:,2)==0 & Orders(:,3)==2);
    
    P_100 = Coeff(Orders(:,1)==1 & Orders(:,2)==0 & Orders(:,3)==0);
    P_010 = Coeff(Orders(:,1)==0 & Orders(:,2)==1 & Orders(:,3)==0);
    P_001 = Coeff(Orders(:,1)==0 & Orders(:,2)==0 & Orders(:,3)==1);
    
    P_110 = Coeff(Orders(:,1)==1 & Orders(:,2)==1 & Orders(:,3)==0);
    P_011 = Coeff(Orders(:,1)==0 & Orders(:,2)==1 & Orders(:,3)==1);
    P_101 = Coeff(Orders(:,1)==1 & Orders(:,2)==0 & Orders(:,3)==1);
    
    J = -1*2*(-4*P_200*P_020*P_002 + ...
           P_200*P_011^2 + P_020*P_101^2 + P_002*P_110^2 + ...
           -1*P_110*P_101*P_011);
       
    if J>=0
        dX = -100;
        dY = -100;
        dZ = -100;    
    else
       
        dX = -4*P_100*P_020*P_002 + ...
              2*P_110*P_010*P_002 + 2*P_101*P_020*P_001 + ...
             -1*P_110*P_011*P_001 - P_101*P_010*P_011 + ...
              P_100*P_011^2;

        dY = -4*P_010*P_200*P_002 + ...
              2*P_011*P_001*P_200 + 2*P_110*P_100*P_002 + ...
             -1*P_101*P_011*P_100 - P_110*P_001*P_101 + ...
              P_010*P_101^2;      

        dZ = -4*P_001*P_200*P_020 + ...
              2*P_010*P_011*P_200 + 2*P_101*P_100*P_020 + ...
             -1*P_110*P_011*P_100 - P_110*P_010*P_101 + ...
              P_001*P_110^2;        
          
        dX = dX / J;
        dY = dY / J;
        dZ = dZ / J;
    end
    
end