function fixedfa =  make_fixedfa_man(sizeXY, Res, Z_arr)
    finalvol_summed = zeros(sizeXY);

    kx = 1:size(finalvol_summed,1);
    ky = 1:size(finalvol_summed,2);

    MultF_X = 1/(length(kx)*Res);
    MultF_Y = 1/(length(ky)*Res);

    CentPos = round((size(finalvol_summed)+1)/2);
    [KX, KY] = ndgrid((kx-CentPos(1))*MultF_X,(ky-CentPos(2))*MultF_Y);
    q2 = KX.^2 + KY.^2 ;
    clear KX KY KZ

    fixedfa_arr = zeros(numel(Z_arr),numel(q2));
    for i = 1:numel(Z_arr)
        fixedfa_arr(i,:) = fatom_vector(sqrt(q2),Z_arr(i));
    end
    fixedfa = mean(fixedfa_arr,1);
end