clc
clear all
close all

range = log10(127);
rangeTick = 127;
dim = 500;
grain = range/dim;
grainTick = rangeTick/dim;

delta = .2;
gamma = 0.3;
beta = 3.5;

i2=-dim:dim;
x2=i2*grain;
x2Tick = i2*grainTick;
p2=delta*gamma+(1-delta)*(1-(1-gamma)*exp(-10.^(beta*x2)));

figure;plot(x2,p2)
set(gca,'xtick',x2(1:100:end),'xticklabel',x2Tick(1:100:end))