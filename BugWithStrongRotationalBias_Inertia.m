clear all;
close all; 
clc;

set(groot,'DefaultTextInterpreter' ,'LaTeX');
set(groot,'DefaultAxesTickLabelInterpreter' ,'LaTeX');
set(groot,'DefaultAxesFontName' ,'LaTeX');
set(groot,'DefaultLegendInterpreter' ,'LaTeX');

steps = 50000000;

rng(654)

h = 0.00005;
v = 3.63463309;
omega = 0.07097372;
Dr = 0.00136558;
D = 0.57998460;
runs = 50;

tau = 1.4;
tau_r = 1.4;

phi = zeros(2,steps);
vec = zeros(4,steps);
t = (0:(steps-1)) * h;

for run = 1:runs

    t = (((run-1)*(steps-1)):((steps-1)*run)) * h;

    vec(:,1) = vec(:,end);
    phi(:,1) = phi(:,end);
    
    for i = 2:steps
       R1 = normrnd(0,1,[3, 1]); 
       phi(:,i) = phi(:,i-1) ...
           + [phi(2,i-1) * h; - phi(2,i-1) * h / tau_r + omega * h / tau_r + sqrt(2 * Dr * h) * R1(1) / tau_r];
       vec(:,i) = vec(:,i - 1) ...
           + [vec(3:4,i-1) * h; - 1 / tau * vec(3:4,i-1) * h + v * h / tau * [ cos(phi(1,i-1)); sin(phi(1,i-1))] + sqrt(2 * D * h) * R1(2:end) / tau];
    end
    
    writematrix(num2str([t(1:10000:end)', vec(1,1:10000:end)', vec(2,1:10000:end)'],'%.6f '),'Data.dat',...
             'Delimiter', 'tab', 'WriteMode', 'append');
end

