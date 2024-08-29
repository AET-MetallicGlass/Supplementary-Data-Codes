% Polynomial tracing based on reconstruction of original orientation

%% load data

SearchRadAr = 3;

Recon_filename = 'Factor_MG_multislice.mat';
filename_interpData = 'finalRec_MG_multislice_intp3.mat';

ourputstring_prefix = 'polyn_tracing_MG_multislice_rad%i';

%% load reconstructed 3D volume
Dsetvol = importdata(Recon_filename);
Dsetvol = Dsetvol.final_Rec;

xx = (1:size(Dsetvol,1)) - round((size(Dsetvol,1)+1)/2);
yy = (1:size(Dsetvol,2)) - round((size(Dsetvol,2)+1)/2);
zz = (1:size(Dsetvol,3)) - round((size(Dsetvol,3)+1)/2);
xxi = ((3*xx(1)):(xx(end)*3))/3;
yyi = ((3*yy(1)):(yy(end)*3))/3;
zzi = ((3*zz(1)):(zz(end)*3))/3;
xxi = xxi(3:end); yyi = yyi(3:end); zzi = zzi(3:end);
[Y,X,Z] = meshgrid(yy,xx,zz);
[Yi,Xi,Zi] = meshgrid(yyi,xxi,zzi);
Dsetvol = interp3(Y,X,Z,Dsetvol,Yi,Xi,Zi,'spline',0);
clear Xi Yi Zi

FinalVol = My_paddzero(Dsetvol,size(Dsetvol)+20);
FinalVol_single = single(FinalVol);
save(filename_interpData,'FinalVol_single','-v7.3');

[ii,jj,kk] = meshgrid(0:4,0:4,0:4);
fitCoeff = [ii(:),jj(:),kk(:),0*kk(:)];
fitCoeff(sum(fitCoeff,2)>4,:) = [];
fitCoeff(max(fitCoeff,[],2) == 4,4) = -1;

for kkkkk = 1:length(SearchRadAr)
    
Th           = 0;
Res          = 0.347/3;
minDist      = 1.75 / Res;
MaxIter      = 14;
CritIter     = 7;
saveInterval = 5000;

SearchRad    = SearchRadAr(kkkkk);
ourputstring = sprintf(ourputstring_prefix,SearchRad);
disp(ourputstring)

%%
se           = strel3d(3);
dilatedBW    = imdilate(FinalVol, se);
maxPos       = find( FinalVol == dilatedBW & FinalVol > Th );
maxVals      = FinalVol(maxPos);
[~,sortInd]  = sort(maxVals,'descend');
maxNum       = min( 100000, length(sortInd) );
maxPos       = maxPos(sortInd(1:maxNum));

fprintf(1,'numpeak = %d \n',length(maxPos));

maxXYZ = zeros(length(maxPos),3);
for i=1:length(maxPos)
    [xx,yy,zz] = ind2sub(size(FinalVol),maxPos(i));
    maxXYZ(i,:) = [xx yy zz];  
end

clear Dsetvol dilatedBW
%%
Q = 0.5;
cropHalfSize = SearchRad;
[X,Y,Z] = ndgrid(-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize);
SphereInd = find( X.^2+Y.^2+Z.^2 <= (SearchRad+0.5)^2 );
XYZdata.X = X(SphereInd);
XYZdata.Y = Y(SphereInd);
XYZdata.Z = Z(SphereInd);

Orders = fitCoeff(:,1:3);
PosArr = zeros(size(maxXYZ));
TotPosArr = zeros(size(maxXYZ));

exitFlagArr = zeros(1, size(maxXYZ,1));

CoeffArr = repmat(fitCoeff(:,4),[1 size(maxXYZ,1)]);

Alpha = 1;

for i=1:size(maxXYZ,1)
    endFlag = 0;
    consecAccum = 0;
    iterNum = 0;
    while ~endFlag    
        iterNum = iterNum + 1;
        if iterNum>MaxIter
          exitFlagArr(i) = -4;
          endFlag = 1;
        end
        cropXind = maxXYZ(i,1) + (-cropHalfSize:cropHalfSize);
        cropYind = maxXYZ(i,2) + (-cropHalfSize:cropHalfSize);
        cropZind = maxXYZ(i,3) + (-cropHalfSize:cropHalfSize);

        cropVol = FinalVol(cropXind,cropYind,cropZind);

        Pos = PosArr(i,:);
        GaussWeight = exp(-1*Alpha*( (X(SphereInd)-Pos(1)).^2 + (Y(SphereInd)-Pos(2)).^2 + (Z(SphereInd)-Pos(3)).^2 ) / cropHalfSize^2 );
        
        fun = @(p,xdata) calculate_3D_polynomial_Rogers(xdata.X,xdata.Y,xdata.Z,Pos,Orders,p).*GaussWeight;

        opts = optimset('Display','off');
        
        [p1,fminres1] = lsqcurvefit(fun,CoeffArr(:,i),XYZdata,cropVol(SphereInd).*GaussWeight,[],[],opts);
        CoeffArr(:,i) = p1;
        
        [dX, dY, dZ] = calc_dX_dY_dZ_Rogers(Orders,CoeffArr(:,i));
        if dX ==-100 && dY == -100 && dZ == -100
            exitFlagArr(i) = -1;
            endFlag = 1;
        else
            maxedShift = max([dX dY dZ],-1*[Q Q Q]);
            minedShift = min(maxedShift,[Q Q Q]);
            PosArr(i,:) = PosArr(i,:) + minedShift;
            if max(abs(PosArr(i,:))) > cropHalfSize
                exitFlagArr(i) = -2;
                endFlag = 1;
            elseif max(abs(minedShift)) < Q
                if consecAccum == CritIter-1
                    goodAtomTotPos = TotPosArr(1:i-1,:);
                    goodAtomTotPos = goodAtomTotPos(exitFlagArr(1:i-1)==0,:);
                    Dist = sqrt(sum((goodAtomTotPos - repmat(PosArr(i,:)+maxXYZ(i,:),[size(goodAtomTotPos,1) 1])).^2,2));
                    if min(Dist) < minDist
                      exitFlagArr(i) = -3;
                    else
                      TotPosArr(i,:) = PosArr(i,:) + maxXYZ(i,:);
                    end
                    endFlag = 1;
                else
                    consecAccum = consecAccum + 1;
                end
            else
                consecAccum = 0;
            end
        end
    end
    fprintf(1,'peak %d, flag %d \n',i,exitFlagArr(i));
    
    if mod(i,saveInterval) == 0
        save(sprintf('%s_result.mat',ourputstring),'PosArr','TotPosArr','CoeffArr','Orders','exitFlagArr');
    end
end

%%
save(sprintf('%s_result.mat',ourputstring),'PosArr','TotPosArr','CoeffArr','Orders','exitFlagArr');

end