% Polynomial tracing based on reconstruction of original orientation
clear all; close all
%% load data
addpath('./src');
homeDir = './';
inputDir = [homeDir, './'];

Recon_filename = [homeDir 'Resire_CoPdPt_Dose1p7e4_rec.mat'] % rec_df_Ge_220326_t1_ref_wonp.mat %ReconObj_Ge_220326_Run1.mat% ReconObj_Ge_220326_Run4000
ourputstring_prefix = [inputDir 'polyn_tracing_CoPdPt_Dose1p7e4']
%% load reconstructed 3D volume
Dsetvol0 = importdata(Recon_filename);
Dsetvol0 = Dsetvol0.reconstruction;
Mask=importdata('./Resire_CoPdPt_Mask.mat');

for i=1:346
    Dsetvol(:,:,i)=Dsetvol0(:,:,i).*Mask(:,:,i);
end

N = size(Dsetvol)
fprintf('Recstruction Dimension = [%d,%d,%d];\n',N(1),N(2),N(3));
%return
Dsetvol0 = My_stripzero(Dsetvol,[N(1),N(2),N(3)]);
Dsetvol_1 = Dsetvol0(:,1:end,:);

FinalVol_tracing = (Dsetvol_1);
N1=size(FinalVol_tracing);

fprintf('Tracing Dimension = [%d,%d,%d];\n',N1(1),N1(2),N1(3));
%% oversamplling the reconstruction matrix by 3*3*3 by linear interpolation
ov=3;% oversamplling
xx = (1:size(FinalVol_tracing,1)) - round((size(FinalVol_tracing,1)+1)/2);
yy = (1:size(FinalVol_tracing,2)) - round((size(FinalVol_tracing,2)+1)/2);
zz = (1:size(FinalVol_tracing,3)) - round((size(FinalVol_tracing,3)+1)/2);

xxi = ((ov*xx(1)):(xx(end)*ov))/ov;
yyi = ((ov*yy(1)):(yy(end)*ov))/ov;
zzi = ((ov*zz(1)):(zz(end)*ov))/ov;

xxi = xxi(3:end); yyi = yyi(3:end); zzi = zzi(3:end);

[Y,X,Z]     = meshgrid(yy,xx,zz);
[Yi,Xi,Zi]  = meshgrid(yyi,xxi,zzi);

FinalVol_tracing=single(FinalVol_tracing);
FinalVol_tracing_interp     = interp3(Y,X,Z,FinalVol_tracing,Yi,Xi,Zi,'spline',0);
FinalVol    = My_paddzero(FinalVol_tracing_interp,size(FinalVol_tracing_interp)+20);% FinalVol_tracing_interp; % My_paddzero(FinalVol_tracing_interp,size(FinalVol_tracing_interp)+20);

FinalVol=double(FinalVol);
N2=size(FinalVol);
fprintf('Tracing Dimension (oversamplling) = [%d,%d,%d];\n',N2(1),N2(2),N2(3));

%%  get polynomial power array
[ii,jj,kk] = meshgrid(0:4,0:4,0:4);
fitCoeff = [ii(:),jj(:),kk(:),0*kk(:)];
fitCoeff(sum(fitCoeff,2)>4,:) = [];
fitCoeff(max(fitCoeff,[],2) == 4,4) = -1;
% return
%%
parpool_size=12;
if parpool_size~=0
    pjob = gcp('nocreate');
    if isempty(pjob)
        parpool(parpool_size)
    elseif pjob.NumWorkers ~= parpool_size
        delete(pjob)
        parpool(parpool_size)
    end
end
%% Tracing parameters

MaxIter = 14;       CritIter = 7;         Th        = 0; 
Res     = 0.347/ov;  minDist  = 1.8/Res;     SearchRad = ov;

%% get the local maxima from the reconstruction volume    
    se           = strel3d(3)
   dilatedBW    = imdilate(FinalVol, se);
    maxPos       = find( FinalVol == dilatedBW & FinalVol > Th );
    maxVals      = FinalVol(maxPos);
    [~,sortInd]  = sort(maxVals,'descend');

    maxNum      = min(100000,length(sortInd));
    maxPos       = maxPos(sortInd(1:maxNum));
    
    fprintf(1,'numpeak = %d \n',length(maxPos));
 
    maxXYZ = zeros(length(maxPos),3);
    for i=1:length(maxPos)
        [xx,yy,zz] = ind2sub(size(FinalVol),maxPos(i));
        maxXYZ(i,:) = [xx yy zz];
    end

    clear Xi Yi Zi Dsetvol dilatedBW
    %%
    Q = 0.5;  Alpha = 1;
    cropHalfSize = SearchRad;
    [X,Y,Z] = ndgrid(-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize);

    SphereInd = find( X.^2+Y.^2+Z.^2 <= (SearchRad+0.5)^2 );
    X_sph = X(SphereInd);
    Y_sph = Y(SphereInd);
    Z_sph = Z(SphereInd);
    
    XYZdata.X = X_sph;
    XYZdata.Y = Y_sph;
    XYZdata.Z = Z_sph;
  
    Orders = fitCoeff(:,1:3);
    PosArr = zeros(size(maxXYZ));
    TotPosArr = zeros(size(maxXYZ));
    
    exitFlagArr = zeros(1, size(maxXYZ,1));
    
    CoeffArr = repmat(fitCoeff(:,4),[1 size(maxXYZ,1)]);

tic
    parfor i=1:size(maxXYZ,1)
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
            GaussWeight = exp(-1*Alpha*( (X_sph-Pos(1)).^2 + (Y_sph-Pos(2)).^2 + (Z_sph-Pos(3)).^2 ) / cropHalfSize^2 ); 
            fun = @(p,xdata) calculate_3D_polynomial_Rogers(xdata.X,xdata.Y,xdata.Z,Pos,Orders,p).*GaussWeight; 
            opts = optimset('Display','off');
            [p1,fminres1] = lsqcurvefit(fun,CoeffArr(:,i),XYZdata,cropVol(SphereInd).*GaussWeight,[],[],opts);
            CoeffArr(:,i) = p1;
%%
            [dX, dY, dZ] = calc_dX_dY_dZ_Rogers(Orders,CoeffArr(:,i));

            if dX ==-100 && dY == -100 && dZ == -100
                exitFlagArr(i) = -1;
                endFlag = 1;
            else
                confinedShift = min(max([dX dY dZ],-1*[Q Q Q]),[Q Q Q]);
                PosArr(i,:) = PosArr(i,:) + confinedShift;
                if max(abs(PosArr(i,:))) > cropHalfSize
                    exitFlagArr(i) = -2;
                    endFlag = 1;
                elseif max(abs(confinedShift)) < Q
                    if consecAccum == CritIter-1  %|| max(abs(confinedShift)) < 1e-5
                        TotPosArr(i,:) = PosArr(i,:) + maxXYZ(i,:); 
                        endFlag = 1;
                    else
                        consecAccum = consecAccum + 1;
                    end
                else
                    consecAccum = 0;
                end
            end
        end
        fprintf(1,'peak %d, flag %d, consecAccum %d\n',i,exitFlagArr(i),consecAccum);
   %              toc
    end

 toc    
%% apply cutoff according to intensity
    for i = 1:size(maxXYZ,1)
        if exitFlagArr(i) == 0
            goodAtomTotPos = TotPosArr(1:i-1,:);
            goodAtomTotPos = goodAtomTotPos(exitFlagArr(1:i-1)==0,:);
            Dist = pdist2(PosArr(i,:)+maxXYZ(i,:),goodAtomTotPos);
            if min(Dist) < minDist
                exitFlagArr(i) = -3;
            end
        end
    end
    save(sprintf('%s_result.mat',ourputstring_prefix),'PosArr','TotPosArr','CoeffArr','Orders','exitFlagArr','SearchRad','maxXYZ');
    


