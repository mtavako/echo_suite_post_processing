function [nBeats,sampInds] = determine_number_of_beats_v2(lvLength,dt,HR)
% Create a vector for data interpolation when used with smoothing
xInit   =   linspace(1,numel(lvLength),numel(lvLength));
% Run spline smoothing on the lvLength waveform
LV_lngth =   fnval(spaps(xInit,lvLength,1E-10),xInit);
% sp = spaps(xInit,lvLength,1E-10);
% LV_lngth2 =   fnval(spaps(xInit,lvLength*100,1E0),xInit);
% figure()
% hold on
% fnplt(sp)
% plot(lvLength)
% plot(LV_lngth)
% plot(LV_lngth2)
% legend('lvLength','LV lngth');
% Find in peaks from mitral annular velocity waveform after spline fitting
[pks,locs]  =   findpeaks(LV_lngth);
% Remove locations where the pks are less than the standard deviation
% around the mean (assuming mean is close to zero), positive peaks only
% LV length is always positive
% !!! mean is not close to zero 
locs(pks < std(LV_lngth(LV_lngth>0))) = [];
% locs(pks < (mean(LV_lngth)+std(LV_lngth(LV_lngth>0)))) = [];
% Remove peaks where the pks are less than the standard deviation
% around the mean (assuming mean is close to zero), positive peaks only
pks(pks  < std(LV_lngth(LV_lngth>0))) = [];
% pks(pks < (mean(LV_lngth)+std(LV_lngth(LV_lngth>0)))) = [];
% Find the distance between peaks
dlocs = diff(locs);
% Set a rejection threshold based on the patient heart rate if it exists;
% if it doesn't, default to 1/10th the length of the signal
if ~isempty(HR)
    thrsh = round(0.25*(60/HR/dt));
else
    thrsh = round(0.1*numel(xInit));
end
% Reject peaks below the threshold if they exist, else keep the signal as
% is and set the peaks to 1 (will indicate only a single beat exists)
try
    pks = pks(logical([1 (dlocs > thrsh)]));
    [pks,~]     =   sort(pks,'descend');
catch
    pks = 1;
end
% Using the located peaks, try to determine the number of signal periods
if numel(pks) > 1
    % Using the lvLength signals, find the signal period length
    % Normalize the lvLenght signal
    nSig1       =   (lvLength-mean(lvLength))/std(lvLength);
    % Autocorrelate the signal
    corSig1     =   fftshift(xcorr(nSig1,nSig1));
    corrSig     =   corSig1(1:round(numel(corSig1)/2));
    % Use findpeaks to find the average period length
    [pks,locs]  =   findpeaks(corrSig);
    [~,inds]    =   sort(pks,'descend');
    locs        =   locs(inds);
    % Automatically compute the number of beats from total signal length
    nBeats      =   round(numel(lvLength)/locs(1));
    % Automatically compute the indices over which measurements will be
    % computed
    sampInds    =   1:locs:numel(lvLength);
    if sampInds(end) < numel(lvLength)
        if   numel(lvLength)-sampInds(end) > 8
            sampInds    =   cat(2,sampInds,numel(lvLength));
        end
    end
else
    nBeats      =   1;
    sampInds    =   ceil(linspace(1,numel(lvLength),nBeats+1));
end