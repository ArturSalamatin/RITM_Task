clc
close all
clear all

fntSize = 14;
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', fntSize);
set(0,'defaultTextFontName', 'Times New Roman')
set(0,'defaultTextFontSize', fntSize)

C = 1; % F
L = 1; % H

U0 = 1; % V
I0 = 0; % A
I_max = 1; % A
q0 = U0/C; % C

omega = sqrt(C/L);

f_path = 'Ritm_Task/';
data = readmatrix([f_path, 'output2.txt']);

t = data(:,1);
q_calc = data(:, 2);
I_calc = data(:, 3);
q_analit = data(:, 5);
I_analit = data(:, 6);

figure(100)
hold on 
box on
grid on
axis([0 30 -1 1])
plot(t*omega, q_calc/q0, 'k--', 'Linewidth', 1)
plot(t*omega, I_calc/I_max, 'b--', 'Linewidth', 1)

plot(t*omega, q_analit/q0, 'k-', 'Linewidth', 1)
plot(t*omega, I_analit/I_max, 'b-', 'Linewidth', 1)

xlabel(['{\itt\omega}, ' char(8211)])
ylabel(['{\itq}/{\itq}_{max}, {\itI}/{\itI}_{max}, ' char(8211)])

figName = 'figure_1A';
% saveas(gca, [figName], 'eps')
saveas(gca, ['Figs\', figName], 'emf')
% saveas(gca, [figName], 'fig')
% saveas(gca, [figName], 'jpg')

data = readmatrix([f_path, 'output4.txt']);

t = data(:,1);
q_calc = data(:, 2);
I_calc = data(:, 3);
q_analit = data(:, 5);
I_analit = data(:, 6);

figure(200)
hold on 
box on
grid on
axis([70 100 -1 1])
plot(t*omega, q_calc/q0, 'k--', 'Linewidth', 1)
plot(t*omega, I_calc/I_max, 'b--', 'Linewidth', 1)

plot(t*omega, q_analit/q0, 'k-', 'Linewidth', 1)
plot(t*omega, I_analit/I_max, 'b-', 'Linewidth', 1)

xlabel(['{\itt\omega}, ' char(8211)])
ylabel(['{\itq}/{\itq}_{max}, {\itI}/{\itI}_{max}, ' char(8211)])

figName = 'figure_1B';
% saveas(gca, [figName], 'eps')
saveas(gca, ['Figs\', figName], 'emf')
% saveas(gca, [figName], 'fig')
% saveas(gca, [figName], 'jpg')

