
function [snd] = generate_soundsRFPc(snd_order, duration, Freq, fade, betw_pause, Fs)
% nur zwei sounds, sin und tria
% order: 1 = sin zuerst, 2 = tria zuerst 
% s_length = länge des sound samples, pro wellenform, in sec
% freq: zB 2000, 2200, 2400
% fade-in und fade-out, in sec
% pause: pause zwischen sounds, in sec
%Fs       = 44000;
% see also https://en.wikipedia.org/wiki/Sawtooth_wave
% 
% To make a sawtooth wave, start with a sine wave,
% 
% sin(2πft) ,
% 
% in which  f  is the frequency and  t  is time. Add to it
% 
% −12sin(2π2ft) 
% 
% (note the minus sign), and
% 
% +13sin(2π3ft) ,
% 
% and so on. The terms have consecutive multiples of the frequency, are divided by that multiple, and have alternating signs.
% 
% Add on an infinite number of terms following that pattern, and you have a sawtooth wave.
% 
% Omit the even-numbered terms, the ones with the minus signs, and you have a square wave.
% 
% Replace the sines in the square wave with cosines, and you have a triangle wave.
% 
% https://www.quora.com/How-does-a-sine-wave-lead-to-a-square-wave-and-a-sawtooth-wave

t         = (0:(1/Fs):duration-(1/Fs))';
fade_samp = length((0:(1/Fs):fade-(1/Fs)));
fadefunc  = ones(length(t), 1);
fadefunc(1:fade_samp) = linspace(0, 1, fade_samp);
fadefunc(end-fade_samp+1:end) = linspace(1, 0, fade_samp);

sinwave = sin(2*pi*Freq*t) .* fadefunc;
sawwave = sawtooth(2*pi*Freq*t, 1) .* fadefunc;
paus    = zeros(Fs * betw_pause, 1);
if snd_order == 1
    snd = [sinwave; paus; sawwave];
elseif snd_order == 2
    snd = [sawwave; paus; sinwave];
end

%theta = linspace(-pi, pi, 360);
fade       = 0.02; %fade-in, fade-out
Freq       = 2200; %in Hz
duration   = 0.8;
betw_pause = 0.4;
snd_order  = 1;
Fs         = 48000; % sampling freq
nodd = 0;
phi=0

t         = (0:(1/Fs):duration-(1/Fs))';
%theta     = linspace(-pi, pi, 360);
wave      = zeros(size(theta));

  for k = 0:nodd
  tmp = (-1)^k * (sin((2*k+1)*pi*Freq*t + phi) / (2*k+1)^2);
  wave = wave + tmp;
  mt(k+1,:)=wave;
  end
sound(wave, Fs)
end

