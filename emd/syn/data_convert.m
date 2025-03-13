clear; close all; clc;

[Data, SegyTraceHeaders, SegyHeader] = ReadSegy('nrand.sgy');

lf = 5;
hf = 120;
dt = 0.004;
N = 3;
verb = 1;

nrand_emd_t = fxemd(Data, lf, hf, dt, N, verb);

figure;

dx = 0.02;
[nz, nx] = size(Data);
taxis = (0: nz) * dt;
xaxis = (0: nx) * dx;

subplot(1, 2, 1);
imagesc(xaxis, taxis, Data);
axis image;
title('Original Data');

subplot(1, 2, 2);
imagesc(xaxis, taxis, nrand_emd_t);
axis image;
title('Processed Data');

WriteSegyStructure('nrand_emd_t.sgy',SegyHeader,SegyTraceHeaders,nrand_emd_t);
