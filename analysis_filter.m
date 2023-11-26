%%%%%%%%%%%%%%%%%%%%%%%% analysis %%%%%%%%%%%%%%%%%
[x,Fs] = audioread("speech8khz.wav");
Hd = equiripple_low_pass().numerator;
freqz(Hd);
saveas(gcf,"Filter_response.jpg");
H0 = Hd(1:2:end); %first polyphase
H1 = Hd(2:2:end); %second polyphase
x0 = x(2:2:end); %if x start at -1, the second sample, is first of xd(2n)
x1 = x(1:2:end); %if x start at -1, the first sample is first of xd(2n-1)
t0 = filter(H0,1,x0);
t1 = filter(H1,1,x1);
vd0 = t0+t1; % multiplying with W
vd1 = t0-t1; % multiplying with W

%%%%%%%%%%%%%%%%%%%%%% synthesis %%%%%%%%%%%%%%%%%%%%%
pre_k0 = vd0+vd1;
pre_k1 = vd0-vd1;
k0 = H1;
k1 = H0;
y0 = filter(k0,1,pre_k0);
y1 = filter(k1,1,pre_k1);

y_good = zeros(length(y0)+length(y1),1);
for i = 1:length(y0)
    y_good(2*i) = y0(i);
end
for i=1:length(y1)
    y_good(2*i-1) = y1(i); %start of commutator
end
temph = zeros(1,length(Hd));
for i=1:length(Hd)
    temph(i) = (-1)^i*(Hd(i));
end
%T = 1/2*(conv(Hd,Hd)-conv(temph,temph));
T = conv(upsample(H0,2),upsample(H1,2));
T = conv(T,[0,1]);
T = 2*T; % T(Z) = 2*Z^(-1)*H0(Z2)*H1(Z2)
zerophase(T,1);
% hold on
% [H,w] = zerophase(Hd);
% plot(w/pi, H.^2);
% plot(1-w/pi, H.^2);
saveas(gcf,"zerophase.jpg");

k0 = H0;
k1 = H1;
y0 = filter(k0,1,pre_k0);
y1 = filter(k1,1,pre_k1);

y_bad = zeros(length(y0)+length(y1),1);
for i = 1:length(y0)
    y_bad(2*i) = y0(i);
end
for i=1:length(y1)
    y_bad(2*i-1) = y1(i);
end

Y_good = fft(y_good);
Y_bad = fft(y_bad);

L = length(y_good);
plot(Fs/L*(-L/2:L/2-1),abs(fftshift(Y_good)),"LineWidth",3)
title("Magnitude spectrum")
xlabel("f(Hz)")
ylabel("|fft(X)|")
saveas(gcf,"y_good.jpg")
L = length(y_bad);
plot(Fs/L*(-L/2:L/2-1),abs(fftshift(Y_bad)),"LineWidth",3)
title("Complex Magnitude of fft Spectrum")
xlabel("f (Hz)")
ylabel("|fft(X)|")
saveas(gcf,"y_bad.jpg")
