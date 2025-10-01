clear all;
close all; 
clc;

set(groot,'DefaultTextInterpreter' ,'LaTeX');
set(groot,'DefaultAxesTickLabelInterpreter' ,'LaTeX');
set(groot,'DefaultAxesFontName' ,'LaTeX');
set(groot,'DefaultLegendInterpreter' ,'LaTeX');

steps = 50000000;

rng(10)

h = 0.0001;
v = 3.63591097;
omega = 0.07095621;
Dr = 0.00124054;
D = 0.65054314;
runs = 500;

phi = zeros(1,steps);
vec = zeros(2,steps);
t = (0:(steps-1)) * h;

for run = 1:runs

    t = (((run-1)*(steps-1)):((steps-1)*run)) * h;

    vec(:,1) = vec(:,end);
    phi(1) = phi(end);
    
    for i = 2:steps
       R1 = normrnd(0,1,[3, 1]); 
       phi(i) = phi(i-1) + omega * h + sqrt(2 * Dr * h) * R1(1);
       vec(:,i) = vec(:,i - 1) + v * [cos(phi(i-1)); sin(phi(i-1))] * h + sqrt(2 * D * h) * R1(2:end);
    end
    
    writematrix(num2str([t(1:10000:end)', vec(1,1:10000:end)', vec(2,1:10000:end)'],'%.6f '),'Data.dat',...
             'Delimiter', 'tab', 'WriteMode', 'append');
end
