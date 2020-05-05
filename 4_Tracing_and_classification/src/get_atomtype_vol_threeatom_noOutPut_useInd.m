function [atomtype] = get_atomtype_vol_threeatom_noOutPut_useInd(...
    DataMatrix, atom_pos, avatom1, avatom2, avatom3, boxhalfsize,useInd)

atomtype = zeros(1,size(atom_pos,2));

for kkk=1:size(atom_pos,2)
    curr_x = round(atom_pos(1,kkk));
    curr_y = round(atom_pos(2,kkk));
    curr_z = round(atom_pos(3,kkk));
    
    curr_vol = DataMatrix(curr_x-boxhalfsize:curr_x+boxhalfsize,...
        curr_y-boxhalfsize:curr_y+boxhalfsize,...
        curr_z-boxhalfsize:curr_z+boxhalfsize);
    
    D_type1 = sum(abs(curr_vol(useInd)-avatom1(useInd)));   % Rfactor N=1  ; all elements > 0
    D_type2 = sum(abs(curr_vol(useInd)-avatom2(useInd)));
    D_type3 = sum(abs(curr_vol(useInd)-avatom3(useInd)));

    D_ar = [D_type1; D_type2; D_type3];
    
    [~, MinInd] = min(D_ar);
    
    if MinInd == 1
        atomtype(kkk) = 1;
    elseif MinInd == 2
        atomtype(kkk) = 2;
    elseif MinInd == 3
        atomtype(kkk) = 3;
    end
end
end