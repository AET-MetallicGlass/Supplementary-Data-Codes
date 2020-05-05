function Q_l_NNN= obtain_Q_l_2ndNN_inds_PBC(Model,BondTh,l,CellPara, PBCEnforceYN)
%% Code for calculating bond length order parameter
% Model: xyz coordinate
% For Pt MD simulation, take BondTh = 3.8
% CellPara, PBCEnforceYN: for periodic bondary condition
%               CellPara: the size of simulation cell
%               PBCEnforceYN: The flag open PBC

Q_l_NNN = zeros(1,size(Model,2));
Q_lm_ar = zeros(size(Model,2),l*2+1);

for i=1:size(Model,2)
    
    currPos = Model(:,i);
    
    SubDist = Model - repmat(currPos,[1 size(Model,2)]);
    if PBCEnforceYN
        SubDist(1,:) = mod(SubDist(1,:) + CellPara(1)/2, CellPara(1)) - CellPara(1)/2;
        SubDist(2,:) = mod(SubDist(2,:) + CellPara(2)/2, CellPara(2)) - CellPara(2)/2;
        SubDist(3,:) = mod(SubDist(3,:) + CellPara(3)/2, CellPara(3)) - CellPara(3)/2;
    end
    
    Dist = sqrt( sum( (SubDist).^2, 1) );
    BondAtomInds = find(Dist < BondTh & Dist > 1e-6);
    
    Q_lm = zeros(length(BondAtomInds),2*l+1);
    
    for j=1:length(BondAtomInds)
        
        currBondVec = currPos - Model(:,BondAtomInds(j));
        
        if PBCEnforceYN
            currBondVec(1,:) = mod(currBondVec(1,:) + CellPara(1)/2, CellPara(1)) - CellPara(1)/2;
            currBondVec(2,:) = mod(currBondVec(2,:) + CellPara(2)/2, CellPara(2)) - CellPara(2)/2;
            currBondVec(3,:) = mod(currBondVec(3,:) + CellPara(3)/2, CellPara(3)) - CellPara(3)/2;
        end

        currBondVec = currBondVec / norm(currBondVec);

        Theta = acos(dot(currBondVec, [0; 0; 1]));
        AsVec = currBondVec - dot(currBondVec, [0; 0; 1])*[0; 0; 1];
        Phi = atan2(AsVec(2),AsVec(1));
        P = legendre(l,cos(Theta));
        P_neg = zeros(l-1,1);
        for N = 1:l
            P_neg(N) = P(l-N+1+1) * (-1)^(l-N+1) * factorial(l-(l-N+1)) / factorial(l+l-N+1);
        end
        
%         P_neg_new = [];
%         for N = 2:l+1
%             P_neg_new(N-1) = P(N) * (-1)^(N-1) * factorial(l-(N-1)) / factorial(l+N-1);
%         end
        
        P = cat(1,P_neg,P);
            
        for m = -1*l:l
            
            Q_lm(j,m+l+1) = (-1)^(mod(m,2)) * sqrt( (2*l+1) / 4 /pi * factorial(l-m) / factorial(l+m)) * P(m+l+1)*exp(1i*m*Phi);
            %Q_lm(j,m+l+1) =  sqrt( (2*l+1) / 4 /pi * factorial(l-m) / factorial(l+m)) * P(m+l+1)*exp(1i*m*Phi);
            
        end
    end
        
    Q_lm_bar = sum(Q_lm,1) / length(BondAtomInds);
        
    Q_lm_ar(i,:) = Q_lm_bar;
    %Q_l(i) = sqrt(4*pi / (2*l+1) * sum(abs(Q_lm_bar).^2));

end



for i=1:size(Model,2)
    
    currPos = Model(:,i);
    
    SubDist = Model - repmat(currPos,[1 size(Model,2)]);
    if PBCEnforceYN
        SubDist(1,:) = mod(SubDist(1,:) + CellPara(1)/2, CellPara(1)) - CellPara(1)/2;
        SubDist(2,:) = mod(SubDist(2,:) + CellPara(2)/2, CellPara(2)) - CellPara(2)/2;
        SubDist(3,:) = mod(SubDist(3,:) + CellPara(3)/2, CellPara(3)) - CellPara(3)/2;
    end
    
    Dist = sqrt( sum( (SubDist).^2, 1) );
    
    BondAtomInds = find(Dist < BondTh);
    
    Q_lm_NNNbond_ar = zeros(length(BondAtomInds),l*2+1);
    
    for j=1:length(BondAtomInds)
        
        Q_lm_NNNbond_ar(j,:) = Q_lm_ar(BondAtomInds(j),:);
    end
        
    Q_lm_bar = sum(Q_lm_NNNbond_ar,1) / size(Q_lm_NNNbond_ar,1);
        
    
    Q_l_NNN(i) = sqrt(4*pi / (2*l+1) * sum(abs(Q_lm_bar).^2));

end




end