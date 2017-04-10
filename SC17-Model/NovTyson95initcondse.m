function y0 = NovTyson95initcondse

nCdc13free = 0.01;
nMPFphos0 = 0.01;
nUbE = 0.0021;
nMPFinactive = 0.01;
nMPFphos2 = 0.307;
%nCdc2free = 0.622;
nCdc25active = 0.026;
nWee1phosl = 0.2;
nWee1active = 0.5;
nWee1phosr = 0.2;
nW = 1.0;
nMik1 = 0.066;
nCdc25total = 0.8;
nDNA = 0.0;
nTaphos = 0;
nmass = 0.6;
nIEphos = 0.001;
%nX = 0.505;
nMPFactive = 0.051;


%% Now define the vector of initial values used in main_yeast.m
y0 = [nCdc13free, nMPFphos0, nMPFinactive, nMPFphos2, nMPFactive,...
    nCdc25active, nWee1active, nWee1phosr, nWee1phosl, nW,...
    nMik1, nIEphos, nUbE, nTaphos, nCdc25total, nDNA, nmass];
