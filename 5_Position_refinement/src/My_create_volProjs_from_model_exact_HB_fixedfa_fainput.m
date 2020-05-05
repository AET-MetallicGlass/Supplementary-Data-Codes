% origin of model should be (0,0,0)
% model is 3xN array, with real unit (A, nm, etc.)
% Res is pixel size, same unit with model

function Projs = My_create_volProjs_from_model_exact_HB_fixedfa_fainput(...
    model, atomtype, Heights, Bfactors, volsize, Res, CropHalfWidth, Angles, fixedfa)


  
model = model ./ Res;
FHeights = Heights;
FWidths = Bfactors / pi^2 / Res^2;

if length(volsize) == 3
    x = (1:volsize(1)) - round((volsize(1)+1)/2);
    y = (1:volsize(2)) - round((volsize(2)+1)/2);
    z = (1:volsize(3)) - round((volsize(3)+1)/2);
elseif length(volsize) == 1
    x = (1:volsize(1)) - round((volsize(1)+1)/2);
    y = x;
    z = x;
else
    error('volsize should be either length 3 or length 1!')
end
    
sizeX = [length(x) length(y)];

inInd = find(model(1,:) >= min(x) & model(1,:) <= max(x) & ...
             model(2,:) >= min(y) & model(2,:) <= max(y) & ...
             model(3,:) >= min(z) & model(3,:) <= max(z));
         
calcModel = model(:,inInd);
calcAtomtype = atomtype(:,inInd);




cropIndRef = -CropHalfWidth:CropHalfWidth;
[cropX,cropY] = ndgrid(cropIndRef,cropIndRef);

Projs = zeros(sizeX(1),sizeX(2),size(Angles,1));
for i=1:size(Angles,1)
  RM1 = MatrixQuaternionRot([0 0 1],Angles(i,1));  
  RM2 = MatrixQuaternionRot([0 1 0],Angles(i,2));
  RM3 = MatrixQuaternionRot([1 0 0],Angles(i,3));
      
  calcModel_rot = (RM1*RM2*RM3)'*calcModel; 
  
  finalvol_padded = zeros( [sizeX, length(Heights)]);
  cenPos = round((size(finalvol_padded)+1)/2);
  
    for j=1:size(calcModel_rot,2)

        currPos = calcModel_rot(1:2,j) + cenPos(1:2)';
        currRndPos = round(currPos);

        cropInd1 = cropIndRef + currRndPos(1);
        cropInd2 = cropIndRef + currRndPos(2);

        CropVol = finalvol_padded(cropInd1,cropInd2,calcAtomtype(j));

        diffPos = currPos-currRndPos;        
        diffPosZ = calcModel_rot(3,j)-round(calcModel_rot(3,j));
        
        gaussCalc = FHeights(calcAtomtype(j))*exp( -1*( (cropX-diffPos(1)).^2 + (cropY-diffPos(2)).^2 )/FWidths(calcAtomtype(j)) );
        gaussZcalc = exp(-1*(cropIndRef - diffPosZ).^2 / FWidths(calcAtomtype(j)) );
        
        finalvol_padded(cropInd1,cropInd2,calcAtomtype(j)) = CropVol + gaussCalc*sum(gaussZcalc);
    end
    
    finalvol_summed = zeros( sizeX);
    if nargin < 9
        kx = 1:size(finalvol_summed,1);
        ky = 1:size(finalvol_summed,2);

        MultF_X = 1/(length(kx)*Res);
        MultF_Y = 1/(length(ky)*Res);

        CentPos = round((size(finalvol_summed)+1)/2);
        [KX, KY] = ndgrid((kx-CentPos(1))*MultF_X,(ky-CentPos(2))*MultF_Y);
        q2 = KX.^2 + KY.^2 ;
        clear KX KY KZ

        fa78 = fatom_vector( sqrt(q2),78 );
        fa28 = fatom_vector( sqrt(q2),28 );

        fixedfa = 0.5*(fa78+fa28);
    end
    for j=1:length(Heights)
        %fa = reshape(fatom_vector( sqrt(q2),AtomNumbers(j)),sizeX);
        CVol = My_stripzero(finalvol_padded(:,:,j),sizeX);
        FVol = my_fft(CVol);
        FVol = FVol .*  reshape(fixedfa,sizeX) ;
        finalvol_summed =finalvol_summed+FVol;
    end

    Vol = real(my_ifft(finalvol_summed));
    Projs(:,:,i) = sum(Vol,3);
end
    
    
    
    
    
    

