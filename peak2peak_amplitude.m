function amp = peak2peak_amplitude(data, fs, t0, xs)
data = squeeze(data);
xs = squeeze(xs);
if size(data,1) > 1 && size(data,2) > 1
    error('Data must have one dimension.');
end
if size(xs,1) > 1 && size(xs,2) > 1
    error('Time vector must have one dimension.');
end

% find mep peak
mepwindow = data(round((t0 + 5)*fs/1000):round((t0 + 25)*fs/1000));
L = length(mepwindow);              % Length of signal
t = xs(round((t0 + 5)*fs/1000)) + (0:L-1)/fs;
[pks, plocs] = findpeaks(mepwindow);
%         mep_peak = t(plocs(1)); % instant of mep peak
pp = pks == max(pks);

% find mep valley
[vls, vlocs] = findpeaks(-mepwindow);
%         mep_valley = t(vlocs(1)); % instant of mep peak
vv = vls == max(vls);
% use this when i need the position of peak and peak intensity separatly
% if ~isempty(plocs) && ~isempty(vlocs)
%     pmax = [find(xs == t(plocs(pp))) mepwindow(plocs(pp))];
%     pmin = [find(xs == t(vlocs(vv))) mepwindow(vlocs(vv))];
% else
%     pmax = [];
%     pmin = [];
% end
if ~isempty(plocs) && ~isempty(vlocs)
    pmax = mepwindow(plocs(pp));
    pmin = mepwindow(vlocs(vv));
else
%     use this to fix one error of peak selection in 7_L_90 - dont know why
%     figure; plot(data);
%     [x, y] = getpts(gca);
%     minmax = [x y];
%     pmax = minmax(2,2);
%     pmin = minmax(1,2);
%     ang = (0:45:315);
%     tx = sprintf('O angulo a ser corrigido eh: %.1f \ndo eletrodo: %s \namplitude:%.5f\n', ang(i), type_sim, pmax-pmin);
%     disp(tx);
    pmax = 0;
    pmin = 0;
end

amp = pmax - pmin;