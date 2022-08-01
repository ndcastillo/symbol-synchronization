% PAM and PCM Script
% by David Castillo, Soledad Villegas

% Delete Cache's Data
close all
clear all
clc

% Exersise
% To Voice's Signal Transmission with SNQ(Quatized) = 35db, Calculate:
% a) Transmission Rate (PCM)
% b) Bandwidth (PCM) in Hertz
% c) Transmission Rate in 64-PAM
% d) Bandwidth (PAM) in Hertz

% Construction of the Signal
% Create parameters to overlap of sinusoidal
% Voice Signal 4kHz

% Signal Construction
% Signals that component a voice signal
A = [3 5 2 1.5];                    % Signal's Signal
fm = [4e3 3.3e3 2.4e3 1.5e3];         % Frequency Signal
wm = 2*pi*fm;                   % Frecquency in rad/s
tm = 1./fm;                     % Time Period

factor = 50;                    % Sample Factor
frecuenciaNyquist = 2*fm;       % Nyquist Rate
fs = factor*frecuenciaNyquist;  % Sample Frequency
ts = 1./fs;                     % Sample Period
rangoDinamico=5;                % Dynamic Range

% Level's Compute

snrQuatizationdB = 35;          %***** SNR (QUATIZED - dB) ***** 
snrQuatization = 10^(snrQuatizationdB/10);
L = sqrt(snrQuatization/3);     % Levels

if L==0 || L < 0
    disp('Fuera del Rango Establecido')
elseif L > 0 && L<=1
    L=1;
elseif L>1 && L<=2
    L=2;
elseif L>2 && L<=4
    L=4;
elseif L>4 && L<=8
    L=8;
elseif L>8 && L<=16
    L=16;
elseif L>16 && L<=32
    L=32;
elseif L>32 && L<=64
    L=64;
elseif L>64 && L<=128
    L=128;
elseif L>128 && L<=256
    L=256;
end

n = log(L)/log(2);                  % Bits per Sample
fprintf('Exersise\n')
disp("Levels = "+L)           % Print Levels
disp("n = "+n+" bits/sample")     % Print Bits per Sample

%% Graphic Voice Signal
numMuestras = min(tm)./min(ts);
numPeriodos = 10;
t = 0:ts:tm*numPeriodos;

for i=1:1:length(A) 
    x(i,:) = A(i)*cos(2*pi*fm(i).*t);
end

F = x(1,:);
figureSinosoides = figure('Name','Continue Signals');
for i=1:1:length(A)      
    subplot(length(A),1,i);
    plotsinusoides = plot(t,x(i,:));
    xlabel('t[s]')
    ylabel('Voltage[V]')
    title('Senoidal Signal')
    hold on;
    F = F + x(i,:);
end

% Graphic Voice Signal
figureSignalVoice = figure('Name','Voice Signal');
plotSignalVoice = plot(t, F);
    plotSignalVoice.LineWidth = 1.5;
    plotSignalVoice = 'red';
    xlabel('t[s]')
    ylabel('Voltage[V]')
    title('Voice Signal')

% Signal Voice Configuration
ts = min(ts);
tm = min(tm);
d=tm;
fs = max(fs);
fm = max(fm);

% Normalization
F = F/max(F);
F = F*rangoDinamico;

% Square Signal
squareSignal = zeros(1,numPeriodos);
squareSignal(1:1)=1;
squareSignal = repmat(squareSignal,1,numMuestras);
F(end)=[];
t(end)=[];

% Sampling
Fsample = F.*squareSignal;

%% Computing Parameters

% Transmission Rate 
R_b = n*frecuenciaNyquist(1);
disp("a) Transmission Rate: R_s= "+R_b + " bps")

% BandWidth PCM
Tb_pcm = 1/R_b;
B_pcm = 1/(2*Tb_pcm);
disp("b) BandWidth (PCM): B_pcm= "+B_pcm+" Hz")


% 64-PAM Rate
k = log(L)/log(2);
R_s = R_b/k;
disp("c) 64-PAM Rate: R_s= "+R_s + " baudios")

% Ancho de banda PAM
Tb_pam = 1/R_s;
B_pam = 1/(2*Tb_pam);
disp("d) BandWidth (PAM): B_pam = "+ B_pam + " Hz")


%% PAM
PAM = [];
k=1;
t_pam1 = 0:ts/numPeriodos:(d*numPeriodos-ts/numPeriodos);
t_pam = 0:Tb_pam/(numPeriodos-1):numPeriodos*length(Fsample)*(Tb_pam/(numPeriodos-1));
t_pam(end)=[];

for i=1:1:length(Fsample)
    for j=1:1:numPeriodos    
        PAM(k)= Fsample(i);
        k=k+1;
    end
end

k=1;
% Retention
Fretention=reshape(Fsample,numPeriodos,[]);
FretentionSignal = [];
for i=1:1:length(Fretention)
    for j=1:1:numPeriodos
        FretentionSignal(k) = Fretention(1,i);
        k=k+1;
    end
end


%% PCM

% Creo un Vector con los niveles de cuantificacion
a = rangoDinamico*2/L;
valoresCuatificacion = -5+a/2:a:5-a/2;

% Quantizing
quatizedSignal = FretentionSignal;
vector = FretentionSignal;
for i=1:1:length(FretentionSignal)
    if FretentionSignal(i) >= valoresCuatificacion(end)
        quatizedSignal(i)= valoresCuatificacion(end);
        vector(i) = L-1;
    elseif FretentionSignal(i) <= valoresCuatificacion(1)
        quatizedSignal(i)=valoresCuatificacion(1);
        vector(i) = 0;
    else
        for j=1:1:L
            if (FretentionSignal(i) > valoresCuatificacion(j) && FretentionSignal(i) < valoresCuatificacion(j) + a/2) || (FretentionSignal(i) < valoresCuatificacion(j) && FretentionSignal(i) > valoresCuatificacion(j) - a/2) 
                quatizedSignal(i) = valoresCuatificacion(j);
                vector(i)=j-1;
            end
        end
    end
end

pcm =reshape(vector,numPeriodos,[]);
pcm_r = pcm(1,:);
pcm_r=dec2bin(pcm_r);

trama=[];
s = [];
k=1;
q=1;

numSamplePoints=10;
for i=1:1:length(pcm_r)
     for j=1:1:n
         s(q) = string(pcm_r(i,j));
         q=q+1;
         for d=1:1:numSamplePoints
             trama(k) = string(pcm_r(i,j));
             k=k+1;
         end
     end
 end

tb=0:Tb_pcm/(numSamplePoints-1):n*numSamplePoints*length(pcm_r)*(Tb_pcm/(numSamplePoints-1));
tb(end)=[];


%% Graphics PCM


% Tags for Coded
tagsDec=0:1:L-1;
tagsBin=dec2bin(tagsDec);
tagsBin=string(tagsBin);
tagsBin=num2cell(tagsBin);

% Coded Signal
figure('Name','CODED SIGNAL');
plot(t,F,t,quatizedSignal, 'LineWidth',1.5);
    yticks(valoresCuatificacion)
    yticklabels(tagsBin)
    style = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',style,'fontsize',8)
    set(gca,'XTickLabelMode','auto')
    title('Coded Signal')
    ylabel('Levels of Voltage [V]')
    xlabel('t[s]')
    grid on;
    
% Digital Data
figureDigital = figure('Name','PCM Signal - Digital Data');
plotDigital = plot(tb,trama);
    plotDigital.LineWidth = 1.5;
    plotDigital.Color='#0D00EB';
    yticks([0 1])
    axis([0 tb(1000) -2 2]);
    title('PCM Signal - Digital Data');
    ylabel('Bit Value');
    xlabel('t[s]');
    grid on;
    grid minor;
    
% BW_pcm
Y=fft(trama);
P2 = abs(Y/length(trama));
P1 = P2(1:(length(trama))/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = B_pcm*numPeriodos*(0:(length(trama)/2))/length(trama);

% PCM Signal - Frequency Domain
figureDigitalFrequencyDomain = figure('Name','PCM Signal - Frequency Domain');
plotDigitalFrequency = plot(f,P1);
    axis([0 36e3 0 0.5]);
    title('PCM Signal - Frequency Domain');
    ylabel('Voltage [V]');
    xlabel('frequency[Hz]');
    grid on;
    grid minor;
    
%% Graphics PAM

%PAM signal
figurePAM = figure('Name','PAM SIGNAL');
plotPAM = plot(t_pam, PAM,'LineWidth',1.5);
    xlabel('t[s]')
    ylabel('Voltage[V]')
    title('PAM SIGNAL')

% BW_pam
Y_pam=fft(PAM);
P2_pam = abs(Y_pam/length(PAM));
P1_pam = P2_pam(1:(length(PAM))/2+1);
P1_pam(2:end-1) = 2*P1_pam(2:end-1);
f_pam = B_pam*numPeriodos*(0:(length(PAM)/2))/length(PAM);

figurePamFrequency = figure('Name', 'PAM Signal - Frequency Domain');
plotPamFrequency = plot(f_pam, P1_pam);
axis([0 12e3 0 0.5]);
    title('PAM Signal - Frequency Domain');
    ylabel('Voltage [V]');
    xlabel('frequency[Hz]');
    grid on;
    grid minor;

% % Quantized Signal
% figure('Name','QUANTIZED SIGNAL');
% plot(t,F,t,quatizedSignal, 'LineWidth',1.5);
%     yticks(valoresCuatificacion)
%     style = get(gca,'XTickLabel');  
%     set(gca,'XTickLabel',style,'fontsize',8)
%     set(gca,'XTickLabelMode','auto')
%     title('Quantized Signal')
%     ylabel('Levels of Quatization [V]')
%     xlabel('t[s]')
% 


%% Digital Data Animation
% close all
% f1=figure('Name','Digital Data');
% for i=1:1:length(tb)
%     plot(tb(1:i),trama(1:i), 'LineWidth',1.5);
%     axis([0 tb(1000) -2 2])
%     grid on;
%     grid minor;
%     pause(0.0001)
% end


%numpoints = 10, entonces se desplaza trama en 5 posiciones
%asumiendo un Tb/2, asi:

bits_zero = zeros(1,5);
length(trama)
delay = [bits_zero trama];
trama2 = [trama bits_zero];

mt = trama2.*delay;
tb2=0:Tb_pcm/(numSamplePoints-1):n*numSamplePoints*length(pcm_r)*(Tb_pcm/(numSamplePoints-1));
tb2 = [tb2 (tb2(end)+Tb_pcm/(numSamplePoints-1))];
tb2 = [tb2 (tb2(end)+Tb_pcm/(numSamplePoints-1))];
tb2 = [tb2 (tb2(end)+Tb_pcm/(numSamplePoints-1))];
tb2 = [tb2 (tb2(end)+Tb_pcm/(numSamplePoints-1))];

% mt
Y=fft(mt);
P2 = abs(Y/length(mt));
P1 = P2(1:(length(mt))/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = B_pcm*numPeriodos*(0:(length(mt)/2))/length(mt);

figureDigital = figure('Name','PCM Signal - Digital Data');
plotDigital = plot(tb2,mt);
    plotDigital.LineWidth = 1.5;
    plotDigital.Color='#0D00EB';
    yticks([0 1])
    axis([0 tb(1000) -2 2]);
    title('PCM Signal - Digital Data');
    ylabel('Bit Value');
    xlabel('t[s]');
    grid on;
    grid minor;
    
% PCM Signal - Frequency Domain
figureDigitalFrequencyDomain = figure('Name','Espectro de m(t)');
plotDigitalFrequency = plot(f,P1);
    axis([0 36e3 0 0.5]);
    title('Espectro en Frecuencia de m(t)');
    ylabel('Voltage [V]');
    xlabel('frequency[Hz]');
    grid on;
    grid minor;
%% Filtro Ideal

[strong_component index_sc] = max(P1(2:end))
f_sincronia = round(f(index_sc+1)+(f(index_sc+2)-f(index_sc+1))/2)
fo = find(f > f_sincronia-200,1)
vector_ff = find(f < f_sincronia+200)
ff = vector_ff(end)
isolated = f(fo:ff)
BPF = zeros(1,length(P1))
BPF(fo:ff) = P1(fo:ff)

figure
plot(f,BPF)
title('Aplicacion del filtro pasa banda');
    ylabel('Voltage [V]');
    xlabel('frequency[Hz]');
    grid on;
    grid minor;
%% Demodulacion
B = 24e3;
Tb = 1/(2*B);


numSamplePoints=10;
k=1;q=1;
yt = [];
val = 0;

%% data clock
vector = [ones(1,numSamplePoints/2) zeros(1,(numSamplePoints/2))]';
vector = [vector vector];
data = vector;
for i=1:1:(length(trama)/(10*2)-1)
    data = [data vector];
end

clock = reshape(data,1,[]);

clock(1)=[]
t_clock = 0:Tb/9:n*numSamplePoints*length(pcm_r)*(Tb_pcm/(numSamplePoints-1));
t_clock(end)=[];t_clock(end)=[];

figureDigital = figure('Name','Señal de Reloj en Receptor Sincronizada');
subplot(3,1,1)
plotDigital = plot(t_clock, clock);
    plotDigital.LineWidth = 1.5;
    plotDigital.Color='#0D00EB';
    yticks([0 1])
    axis([0 tb(1000) -2 2]);
    title('Señal de Reloj en Receptor Sincronizada');
    ylabel('Bit Value');
    xlabel('t[s]');
    grid on;
    grid minor;

subplot(3,1,2)
plotDigital = plot(tb, trama);
    plotDigital.LineWidth = 1.5;
    plotDigital.Color='#0D00EB';
    yticks([0 1])
    axis([0 tb(1000) -2 2]);
    title('Señal Recibida m(t)');
    ylabel('Bit Value');
    xlabel('t[s]');
    grid on;
    grid minor;

%% sampling
trama(end)=[]
tb(end)=[]
samples = trama.*clock;

k = 1;
for i=1:1:length(samples)
   if mod(i,10) == 1
       yt(k) = samples(i)
       k = k + 1
   end
end

subplot(3,1,3)
plotDigital = plot(tb, samples);
    plotDigital.LineWidth = 1.5;
    plotDigital.Color='#0D00EB';
    yticks([0 1])
    axis([0 tb(1000) -2 2]);
    title('Muestreo y_i');
    ylabel('Bit Value');
    xlabel('t[s]');
    grid on;
    grid minor;

 % Detector de Simbolo

simbolos = string(reshape(yt,n,[])');
dim_simbolos = size(simbolos);
vector_sim =[]
volts_simbol = []

for i=1:dim_simbolos(1)
    sis=[]
    for j=1:dim_simbolos(2)
        sis = strcat(sis,simbolos(i,j));
    end
    vector_sim(i) = sis;
    volts_simbol(i)=valoresCuatificacion(bin2dec(sis)+1);
end

m_estimada=[]
k=1
for i=1:1:length(volts_simbol)
     for d=1:1:numSamplePoints
         m_estimada(k) = volts_simbol(i);
         k=k+1;
     end
end
 
t = 0:ts:tm*numPeriodos;
t(end)=[]

figure
plot(t,m_estimada, 'LineWidth',1.5);
    yticks(valoresCuatificacion)
    yticklabels(tagsBin)
    style = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',style,'fontsize',8)
    set(gca,'XTickLabelMode','auto')
    title('Señal Estimada')
    ylabel('Simbolos')
    xlabel('t[s]')
    grid on;
    
