clear variables
showBool        = false;
Gain            = 0.33;
ContrastFactor  = 0.95;
LineThresh      = 0.01; % ¡¡ No modificar !!

xPix_Window     = 178:220;
xPix_BackWindow = [xPix_Window(1)-floor(size(xPix_Window,2)/2)  :xPix_Window(1)-1                           ,...
                   xPix_Window(end)+1                           :xPix_Window(end)+ceil(size(xPix_Window,2)/2)];
yPix_Window     = 1:2050;
%_____________________________________________________________________________________
sBias = FITS.read2sim('bias\*.fits');
BiasCube = zeros(numel(yPix_Window),numel(xPix_Window),numel(sBias));
for index = 1:numel(sBias)
    sBias(index).Im     = fitsread(sBias(index).ImageFileName,'IMAGE');
    BiasCube(:,:,index) = sBias(index).Im(yPix_Window,xPix_Window)*Gain;
end
MBias = median(BiasCube,3);
if showBool
    figure
    imshow(imrotate(MBias,90),[3304 3311])
    title('MasterBias')
end
%_____________________________________________________________________________________
sFlat = FITS.read2sim('flat\*.fits');
FlatCube = zeros(numel(yPix_Window),numel(xPix_Window),numel(sFlat));
for index = 1:numel(sFlat)
    sFlat(index).Im     = fitsread(sFlat(index).ImageFileName,'IMAGE');
    FlatCube(:,:,index) = sFlat(index).Im(yPix_Window,xPix_Window)*Gain - MBias;
end
MFlat = median(FlatCube,3);
MNormFlat = MFlat/mean(MFlat(:));
if showBool
    figure
    imshow(MNormFlat,[median(MNormFlat(:))*ContrastFactor median(MNormFlat(:))/ContrastFactor])
    title('MasterFlat')
end
%_____________________________________________________________________________________
sHe         = FITS.read2sim('cal\ALyb100428.fits');
sHe.Im      = fitsread(sHe.ImageFileName,'IMAGE');
sHe.Im      = (sHe.Im(yPix_Window,xPix_Window)*Gain - MBias)./MNormFlat;
HeSpec      = flip(sum(sHe.Im,2));

sNe         = FITS.read2sim('cal\ALyb100429.fits');
sNe.Im      = fitsread(sNe.ImageFileName,'IMAGE');
sNe.Im      = (sNe.Im(yPix_Window,xPix_Window)*Gain - MBias)./MNormFlat;
NeSpec      = flip(sum(sNe.Im,2));

sThAr       = FITS.read2sim('cal\ALyb100430.fits');
sThAr.Im    = fitsread(sThAr.ImageFileName,'IMAGE');
sThAr.Im    = (sThAr.Im(yPix_Window,xPix_Window)*Gain - MBias)./MNormFlat;
ThArSpec    = flip(sum(sThAr.Im,2));

if showBool
    figure
    subplot(311)
    imshow(imrotate(sHe.Im,90),[median(sHe.Im(:))*ContrastFactor 10*median(sHe.Im(:))/ContrastFactor])
    title('He lamp')
    subplot(312)
    imshow(imrotate(sNe.Im,90),[median(sNe.Im(:))*ContrastFactor 10*median(sNe.Im(:))/(ContrastFactor)])
    title('Ne lamp')
    subplot(313)
    imshow(imrotate(sHe.Im,90),[median(sHe.Im(:))*ContrastFactor 10*median(sHe.Im(:))/ContrastFactor])
    title('ThAr lamp')
end
%_____________________________________________________________________________________
% Reconocimiento manual de las lineas de Helio
HeLines = find(islocalmax(HeSpec,'MinProminence',min(HeSpec) + LineThresh*(max(HeSpec) - min(HeSpec)),'FlatSelection','center'));
if showBool
    for i=1:numel(HeLines)
        figure(1)
        clf
        plot(HeSpec)
        hold on
        plot(HeLines(i),HeSpec(HeLines(i)),'rx')
        title('He Lamp')
        axis([HeLines(i)-400 HeLines(i)+400 min(HeSpec) max(HeSpec)])
        keyboard
    end
end
HeRefLinesAngstrom  = [3888.646, 3964.727, 4026.189, 4471.477, 4713.143, 4921.929, 5015.675, 5875.618, 6678.149, 7065.188, 7281.349];
HeRefLinesPix       = HeLines;
% La última linea detectada por el script no tiene correspondencia en la grafica, por lo que se descarta:
HeRefLinesPix(end) = [];

%_____________________________________________________________________________________
% Reconocimiento manual de las lineas de Neon
NeLines = find(islocalmax(NeSpec,'MinProminence',median(NeSpec) + LineThresh*(max(NeSpec) - min(NeSpec)),'FlatSelection','center'));
if showBool
    for i=1:numel(NeLines)
        figure(1)
        clf
        plot(NeSpec)
        hold on
        plot(NeLines(i),NeSpec(NeLines(i)),'rx')
        title('Ne Lamp')
        axis([NeLines(i)-400 NeLines(i)+400 min(NeSpec) max(NeSpec)])
        keyboard
    end
end
NeRefLinesAngstrom = [5852.487, 5881.895, 5944.834, 6029.9971, 6074.337, 6096.163, 6143.0623, 6163.5939, 6217.2813, 6266.495, 6334.4279, 6402.246, 6506.5279, 6532.8824, 6598.9529, ...
    6717.0428,6929.468, 7032.4127, 7173.939, 7245.167, 7438.899, 8377.607, 8495.36];
NeRefLinesPix       = NeLines;
% Las lineas 4,12,14,19,26,29,30 no tienen correspondencia en la grafica o estan superpuestas a otras lineas, por lo que se descartan:
NeRefLinesPix([4,12,14,19,26,29,30]) = [];

%_____________________________________________________________________________________
% Reconocimiento manual de las lineas de Torio/Argon             
ThArLines = find(islocalmax(ThArSpec,'MinProminence',median(ThArSpec) + 10*LineThresh*(max(ThArSpec) - min(ThArSpec)),'FlatSelection','center')); % Umbral inferior modificado debido a la cantidad de maximos locales encontrados
if showBool
    for i=1:numel(ThArLines)
        figure(1)
        clf
        plot(ThArSpec)
        hold on
        plot(ThArLines(i),ThArSpec(ThArLines(i)),'rx')
        title('ThAr Lamp')
        axis([ThArLines(i)-400 ThArLines(i)+400 min(ThArSpec) max(ThArSpec)])
        keyboard
    end
end
ThArRefLinesAngstrom = [6965.4307, 7067.2181, 7383.9805, 7635.106, 7948.1764, 8264.5225, 8521.4422,9122.9674];
ThArRefLinesPix       = ThArLines;
% Las lineas 1,5,7,9,10,12 no tienen correspondencia en la grafica o estan superpuestas a otras lineas, por lo que se descartan:
ThArRefLinesPix([1,5,7,9,10,12]) = [];

if showBool
    figure
    plot(HeRefLinesPix,HeRefLinesAngstrom,'mo')
    hold on
    plot(NeRefLinesPix,NeRefLinesAngstrom,'bo')
    plot(ThArRefLinesPix,ThArRefLinesAngstrom,'go')
    title('Relacion PixRow <=> Angstromg')
    xlabel('PixRow')
    ylabel('Angstrom')
    legend({'He ref','Ne ref','ThAr ref'})
end

PixRef      = [HeRefLinesPix      ; NeRefLinesPix      ; ThArRefLinesPix     ];
LambdaRef   = [HeRefLinesAngstrom , NeRefLinesAngstrom , ThArRefLinesAngstrom]';

f = fit(PixRef,LambdaRef,'poly4');
if showBool
    figure
    plot(f,PixRef,LambdaRef)
    legend({'Referencias','Funcion cubica extrapolada'})
    xlabel('PixRow')
    ylabel('Angstrom')
end

%_____________________________________________________________________________________
sSpec       = FITS.read2sim('im\ALyb100251.fits');
sSpec.Im    = fitsread(sSpec.ImageFileName,'IMAGE');
BackIm      = (sSpec.Im(yPix_Window,xPix_BackWindow)*Gain - MBias)./MNormFlat;
Im          = (sSpec.Im(yPix_Window,xPix_Window)*Gain - MBias)./MNormFlat;

BackCol = median(BackIm,2);
Spec1D  = flip(sum(Im,2) - size(Im,2)*BackCol);

SpecLength = numel(Spec1D);
LambdaMin = f(1);
LambdaMax = f(SpecLength);
LambdaResolution = (LambdaMax - LambdaMin)/SpecLength; % Angstrom / pixel;

xData = (LambdaMin:LambdaResolution:LambdaMax-LambdaResolution);
yData = Spec1D;
if showBool
    figure
    plot(xData,yData)
    grid on
    xlabel('Angstrom')
    ylabel('sum(e^-)')
end

TopThresh = islocalmax(medfilt1(yData));
f=fit(xData(TopThresh)',yData(TopThresh),'smoothingspline','SmoothingParam',1e-5);

if showBool
    figure
    yDataNorm = yData./f(xData);
    xDataNorm = xData;
    plot(xDataNorm,yDataNorm,'k')
    axis([3800 5200 0.5 1.1])
    xlabel('Angstrom')
    ylabel('Flujo normalizado')
    title('Espectro normalizado')
end