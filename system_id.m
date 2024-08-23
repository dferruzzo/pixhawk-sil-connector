%% Identificação do sistema
% Carrega os dados 
% pwm.mat
% rates.mat
%
% Os dados são gerados a partir de um escrip em MAVSDK-Python
% utilizando o Offboard mode.
%
% O script processa dos dados para obter: 
% 1. A resposta em frequencia para p(t)
% 2. Identifica o modelo linear p(s)/tau_p(2).
%
% -------------------------------------------------------
clear
close all
clc
% Carregar dados
pwm = load('pwm.mat');
time = pwm.ans.Time;
motor1 = pwm.ans.Data(:,1);
motor2 = pwm.ans.Data(:,2);
motor3 = pwm.ans.Data(:,3);
motor4 = pwm.ans.Data(:,4);
%
rates = load('rates.mat');
time_pqr = rates.ans.Time;
p = rates.ans.Data(:,1);
q = rates.ans.Data(:,2);
r = rates.ans.Data(:,3);
% Gráficos
figure(1);
plot(time,motor1,...
     time,motor2,...
     time,motor3,...
     time,motor4,...
     time,p,...
     time,q,...
     time,r);
legend('m1','m2','m3','m4','p','q','r');
xlabel('tempo (s)');
grid on;
%% Processamento dos sinais
% Recortar os sinais
ti_clip = 19.0; % segundos
tf_clip = 26.4; % segundos
idx_ti = find(time<=ti_clip);
idx_tf = find(time>=tf_clip);
t_clipped = time(idx_ti(end):idx_tf(1));
m1_clipped = motor1(idx_ti(end):idx_tf(1));
m2_clipped = motor2(idx_ti(end):idx_tf(1));
m3_clipped = motor3(idx_ti(end):idx_tf(1));
m4_clipped = motor4(idx_ti(end):idx_tf(1));
p_clipped = p(idx_ti(end):idx_tf(1));
q_clipped = q(idx_ti(end):idx_tf(1));
r_clipped = r(idx_ti(end):idx_tf(1));
% Tirando a média dos motores
m1_clipped_no_mean = m1_clipped - mean(m1_clipped);
m2_clipped_no_mean = m2_clipped - mean(m2_clipped);
m3_clipped_no_mean = m3_clipped - mean(m3_clipped);
m4_clipped_no_mean = m4_clipped - mean(m4_clipped);
% juntando os sinais dos motores m1 + m2 - m3 - m4
input = m2_clipped_no_mean + ...
        m3_clipped_no_mean - ...
        m1_clipped_no_mean - ...
        m4_clipped_no_mean;
%%
% Gráficos
figure(2);
plot(t_clipped, input,...
     t_clipped, p_clipped);
legend('input','p');
xlabel('tempo (s)');
grid on;
%% FFT dos sinais
% tempo de amostragem
T = time(2)-time(1); 
% frequencia de amostragem
Fs = 1/T;
% tamanho do sinal
L = length(t_clipped);
% single-sided spectrum da entrada
fft_input_2s = abs(fft(input)/L);
fft_input = fft_input_2s(1:L/2+1);
fft_input(2:end-1) = 2*fft_input(2:end-1);
freqs = Fs/L*(0:(L/2));
freqs_rad = freqs*2*pi;
%
% single-sided spectrum de p
fft_p_2s = abs(fft(p_clipped)/L);
fft_p = fft_p_2s(1:L/2+1);
fft_p(2:end-1) = 2*fft_p(2:end-1);
%
figure(3);
plot(freqs_rad, fft_input,...
     freqs_rad, fft_p); 
title("Single-Sided Amplitude Spectrum");
xlabel("w (rad/s)");
legend('fft input','ffp p');
grid on;
%% Diagrama de Bode de amplitude e fase
% Resposta em frequencia de amplitude de p(t) (rfa_p)
rfa_p = NaN*ones(length(fft_p),1);
for i=1:length(fft_p)
    rfa_p(i) = 20*log10(abs(fft_p(i))/abs(fft_input(i)));
end

Yi = fft(input);
Yi = Yi(1:L/2+1);
Yo = fft(p_clipped);
Yo = Yo(1:L/2+1);
%
%tol = 1e-6;
%Yi(abs(Yi) < tol)=NaN;
%Yo(abs(Yo) < tol)=NaN;
%
angle_dif = (unwrap(angle(Yo))-unwrap(angle(Yi)))*180/pi;
%angle_dif = unwrap(angle(Yo));
%% figures Diagramas de Bode de p(t)
figure(4);
subplot(2,1,1);
semilogx(freqs_rad,rfa_p);
title('Resposta em frequencia de amplitude p(t)');
ylabel('20log_{10}|A_o/A_i| (dB)');
xlim([0, 200])
%xlim([0 120]);
%xlabel('w (rad/s)');
grid on;
subplot(2,1,2);
semilogx(freqs_rad, angle_dif);
title('Resposta em frequenciade fase p(t)');
%ylim([-2000,0]);
xlim([0, 200])
xlabel('w (rad/s)');
grid on;
%% Guradando os dados
resp_freq_p_vector = [freqs_rad', rfa_p, angle_dif];  % Because magic is essential
save('resp_freq_p.mat', 'resp_freq_p_vector');



