fade       = 0.02; %fade-in, fade-out
Freq       = 2200; %in Hz
duration   = 0.8;
betw_pause = 0.4;
snd_order  = 1;
Fs         = 48000; % sampling freq
nodd = 1.5;


snd = generate_sounds(snd_order, duration, Freq, fade, betw_pause, Fs);
sound(snd, Fs);

