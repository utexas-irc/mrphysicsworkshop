function mxyVec = profileLoop(B1, t, pulseshape, BW)

fVec = 0:BW/400:BW/2;

for ii = 1:numel(fVec)
[M T] = pulsesim(B1, t, fVec(ii), pulseshape, -1);
theta(ii) = str2double(T);
Mxy(ii) = abs(M(1,end) + 1i*M(2,end));
end

bwVec  = [-fliplr(fVec) fVec(2:end)];
mxyVec = [fliplr(Mxy)  Mxy(2:end)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find FWHM and datapoints to label %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ltmPts = find(Mxy < max(Mxy)/2); %less than max points
%plotPts = [numel(Mxy)-ltmPts(1) numel(Mxy)+ltmPts(1)];
fwhmX = interp1(Mxy(ltmPts(1))-2:Mxy(ltmPts(1))+2, ...
             fVec(ltmPts(1))-2:fVec(ltmPts(1))+2, ...
             0.5, 'spline');

%fwhmTxt = sprintf('FWHM = %0.2f Hz', bwVec(plotPts(2)) - bwVec(plotPts(1)));
fwhmTxt = sprintf('FWHM = %0.2f Hz', 2*abs(fwhmX));

figure 
plot(bwVec, mxyVec);
hold on
plot([-fwhmX fwhmX], [0.5 0.5], '-ro');
%plot(bwVec(plotPts), mxyVec(plotPts), '-ro');
% fwhmLabel = text(bwVec(plotPts(2))+0.1*numel(bwVec), ...
%                  mxyVec(plotPts(2)),   ...
%                  ['\leftarrow' fwhmTxt], 'FontSize', 12);
fwhmLabel = text(fwhmX+0.1*numel(bwVec), ...
                 0.5,   ...
                 ['\leftarrow' fwhmTxt], 'FontSize', 12);

xlabel('Off-Resonant Frequency (Hz)', 'FontSize', 14);
ylabel('M_{xy} (a.u.)', 'FontSize', 14);
title(['Slice Profile of ' num2str(theta(1)) '\circ ' pulseshape ' Pulse'], 'FontSize', 16, 'FontWeight', 'Bold');




end