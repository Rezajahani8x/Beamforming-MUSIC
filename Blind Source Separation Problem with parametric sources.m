M = 10;
N = 2;
d = 1;
fs = 1e6;
t_T = 1e-3;
f1 = 20e3;
f2 = 10e3;
fc = 150e6;
c = 3e8;
k = 2*pi*fc/c;
doa = [10 20];
doa = deg2rad(doa);
t = 0:1/fs:t_T-1/fs;
S = [exp(1i*2*pi*f1*t);exp(1i*2*pi*f2*t)];
W = 1/sqrt(2) * (randn(M,length(t)) + 1i*randn(M,length(t)));
steering_arg = transpose(1i*k*(0:M-1)*d) * sin(doa);
A = exp(steering_arg);
X = A*S;
Y = X + W;

    %% Initial Steps for Q2 and Q3
theta_start = 0;
theta_end = 90;
[U,G,V] = svd(Y);
% V = transpose(V);
D = G(1:M,1:M);
Usig = U(:,1:N);
Unull = U(:,N+1:end);
Vsig = V(1:N,:);
Vnull = V(N+1:end,:);
theta_deg = theta_start:theta_end;
theta = deg2rad(theta_deg);
a = exp(transpose(1i*k*(0:M-1)*d) * sin(theta));

    %% Q2) Beamforming - DOA
spec_beamf = sum(abs((ctranspose(a)*Usig)).^2,2);
% [~,locs] = findpeaks(spec_beamf,'MinPeakHeight',0.5*max(spec_beamf));
[~,locs] = maxk(spec_beamf,N);
estimated_doa_beamf = theta_deg(locs(1:N));
disp("Estimated DOA (Beamforming Method):");
disp(estimated_doa_beamf);

    %% Q3) MUSIC - DOA
spec_music = 1./sum(abs((ctranspose(a)*Unull)).^2,2);
% [~,locs] = findpeaks(spec_music,'MinPeakHeight',0.1*max(spec_music));
[~,locs] = maxk(spec_music,N);
estimated_doa_music = theta_deg(locs(1:N));
disp("Estimated DOA (MUSIC Method):");
disp(estimated_doa_music);

figure(1);
plot(theta_deg,spec_beamf/max(spec_beamf)); hold on; grid on;
plot(theta_deg,spec_music/max(spec_music));
xlabel("\theta");
title("Normalized MUSIC & Beamforming Spectrums for DOA Estimation");
legend("Beamforming","MUSIC");

    %% Initial Steps for Q4 and Q5
f_start = 1;                                        % Per KHz
f_end = 40;                                         % Per KHz
freq = (f_start:f_end);
f = freq * 1e3;
s = exp(1i*2*pi*transpose(f)*t);

    %% Q4) Beamforming - Source Signal
spec_beamf_s = sum(abs(conj(s)*transpose(Vsig)).^2,2);
% [~,locs] = findpeaks(spec_beamf_s,'MinPeakHeight',0.5*max(spec_beamf_s));
[~,locs] = maxk(spec_beamf_s,N);
estimated_f_beamf = f(locs(1:N));
disp("Estimated Frequencies [KHz] (Beamforming Method):");
disp(estimated_f_beamf*1e-3);

    %% Q5) MUSIC - Source Signal
spec_music_s = 1./sum(abs(conj(s)*transpose(Vnull)).^2,2);
% [~,locs] = findpeaks(spec_music_s,'MinPeakHeight',0.8*max(spec_music_s));
[~,locs] = maxk(spec_music_s,N);
estimated_f_music = f(locs(1:N));
disp("Estimated Frequencies [KHz] (MUSIC Method):");
disp(estimated_f_music*1e-3);

figure(2);
plot(freq,spec_beamf_s/max(spec_beamf_s)); hold on; grid on;
plot(freq,(spec_music_s-mean(spec_music_s))/max(spec_music_s));
xlabel("f [KHz]");
title("Normalized MUSIC & Beamforming Spectrums for Frequency of Sources");
legend("Beamforming","MUSIC");