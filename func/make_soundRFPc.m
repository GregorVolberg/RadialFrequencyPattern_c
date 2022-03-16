
function [snd] = make_soundRFPc(duration, Freq, fade, Fs, nodd)
% for ref see https://www.quora.com/How-does-a-sine-wave-lead-to-a-square-wave-and-a-sawtooth-wave
% and https://en.wikipedia.org/wiki/Sawtooth_wave
t         = (0:(1/Fs):duration-(1/Fs))';
fade_samp = length((0:(1/Fs):fade-(1/Fs)));
fadefunc  = ones(length(t), 1);
fadefunc(1:fade_samp) = linspace(0, 1, fade_samp);
fadefunc(end-fade_samp+1:end) = linspace(1, 0, fade_samp);

twav = [];

% subtract even, add uneven
% use values of  1 3 5 7 9 for nodd
for k = 1:nodd
twav(k, :) = sign(mod(k, 2) - 0.5)/k * (sin(2*pi*k*Freq*t));
end
snd = sum(twav, 1);
snd = snd * (max(twav(1,:))/max(snd)); % re-scale
snd = snd .* fadefunc';
snd = repmat(snd,2,1);
end
