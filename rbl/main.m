%%
clear, clc;

%============================= Analisis Neutronik =================================
aq=10000000;
%NEUTRONIK RADIAL
r=170; %Diameter teras (cm)qswa
N=200; %Jumlah Partisi
psirad=zeros(N,1);
rcell=zeros(N,1);
dr=r/N; %jarak antar partisi, asumsikan sama
rcell(1,1)=-r/2;

for i = 2:N
    rcell(i,1)=rcell(i-1,1)+dr;
end

for i = 1:N
    psirad(i,1)=besselj(0,((2.405*rcell(i,1))/(r+30)));
end

sumpsirad=0;
for i = 1:N
    sumpsirad=sumpsirad+psirad(i,1);
end

powcellconsrad=aq/sumpsirad;
for i = 1:N
    powcellrad(i,1)=powcellconsrad*psirad(i,1);
end

%NEUTRONIK AXIAL
D = 15.21; %Konstanta Difusi (cm)
Sa = 183.0681; %0.0924 biasa, jumlah(1.492*densitas)
vSf = 188.2194; %0.095 biasa,jumlah(0.069*densitas)563.5863, v=7%*3.10^8, vsf = 0.483
d = 365.8; %Tinggi teras
dz = d/N; %Jarak antar partisi, asumsikan sama
error = 0.0001; %Error fluks dan Keff
iter = 0; %Inisialisasi

%Dilakukan penebakan secara random dengan menggunakan variabel rand
%Tebak nilai keff dari 1 sampai 10
keff = rand*10;

%Dilakukan pengulangan untuk mengisi tebakan fluks dengan N partisi
%Definisikan matrix fluks psi, ukuran Nx1
for i = 1:N
    psiax(i,1) = rand*10; %Dilakukan penebakan nilai fluks antara 1 sampai 100
end

%Pembuatan matriks A
%Karena partisi sama dan konstanta difusi konstan maka alpha(a) dan beta(b) konstan
A = zeros(N,N); % Inisialisasi matriks A NxN dengan nilai 0 semua
X = zeros(N,1); % Inisialisasi matrix X Nx1 dengan nilai 0 semua
a = -(D/dz); %Alpha
b = 2*(D/dz) + Sa*dz; %Beta

%Pengisian alpha dan beta pada matriks A tridiagonal
for i = 1:N
    if (i == 1)
        A(i,i) = b;
        A(i,i+1) = a;
    elseif (i == N)
        A(i,i) = b; 
        A(i, i-1) = a;
    else
        A(i,i) = b; %%A(i,1) = b;
        A(i, i-1) = a;
        A(i,i+1) = a;
    end
    X(i,1) = i*dz; %Variabel posisi, untuk plot flux terhadap posisi
end

%Iterasi untuk melakukan perhitungan fluks dan Keff
loop = true; %Inisialisasi loop
while(loop)
    %Mengisi matrix source untuk setiap flux
    for i = 1:N
        S(i,1) = (vSf*psiax(i,1)*dz)/keff;
    end
    
    %Menghitung flux baru dari invers A dikali S, matrix lama ditampung di variabel psix
    psiax1 = psiax;
    psiax = A\S;
    
    %Menghitung Fission product F menggunakan flux lama dan baru (Psix dan Psi)
    F = 0; %Fission product baru(Inisialisasi)
    F1 = 0; %Fission product lama (Inisialisasi)
    %Dilakukan penjumlahan untuk setiap fluks
    for i = 1:N
        F1 = F1 + vSf*psiax1(i,1)*dz;
        F = F + vSf*psiax(i,1)*dz;
    end
    
    %Menentukan keff baru, keff lama disimpan ke variabel keff1
    keff1 = keff;
    keff = (F/F1)*keff1;
    
    %Periksa konvergensi dari keff dan Fluks
    ek = abs((keff - keff1)/(keff1));
    ep = abs((psiax - psiax1)\(psiax1));
    
    %Dilakukan pengecekan, jika sudah lebih kecil dari error, looping berhenti
    if((ek < error)&&(ep < error))
        loop = false; %Loop berhenti
    end
    
    %Pengecekan divergen pada keff dan Flux, tebak ulang fluks dan keff
    %Sehingga perhitungan akan berulang dari awal
    if (isnan(ek))
        keff = rand*10;
        for i = 1:N
            psiax(i,1) = rand*10;
        end
    end
    iter = iter+1;
end

%Menampilkan nilai keff 
disp('K_eff : ');
disp(keff);

%{
%======================= Perhitungan Daya =================================
%Inisialisasi
sumpsiax = 0;
powcellax = zeros(N,1);
psiaxcons = zeros(N,1);
tcell = zeros(N,1);

%Perhitungan daya dengan normalisasi
%Normalisasi fluks neutron aksial
for i=1:N
    sumpsiax = sumpsiax+psiax(i,1);
end

%Normalisasi fluks neutron radial
for i=1:N
    psiaxcons(i,1)=powcellrad(i,1)/sumpsiax;
end

%Memperoleh matriks yang berisikan daya per satuan luas
for i=1:N
    for j=1:N
        powcellax(j,i)=psiaxcons(i,1)*psiax(j,1);
    end
end
%}

%{
%======================= Analisis Burn-Up =================================
% Semua nilai cross-section (variabel dengan awalan "scXX/saXX/sfXX") 
% memiliki satuan (atom/barn.cm)
% dan semua half-time (variabel dengan awalan "htXX" memiliki satuan (s)
flux= 0.964613E+13; % neutron/cm^2.s
n   = 0.5;

%detik*menit*jam*hari*bulan*tahun
tf  = 60*60*24*30*12*2; % detik (Waktu = 2 tahun)
dt  = 60*60; % detik
t   = 0:dt:tf;
% t   = [0:dt:tf];

% Fuel
% U-238
Nu8(1)  = 2.2968545203386600E-02;   
scu8    = 2.68e-24;
sfu8    = 16.8e-30;
sau8    = scu8+sfu8;

% U-239
Nu9(1)  = 0;
htu9    = 23.45*60;
lu9     = log(2)/htu9;

% Np-239
Nnp9(1) = 0;
htnp9   = 2.356*24*60*60;
lnp9    = log(2)/htnp9;

% Pu-239
Npu9(1) = 0;
sapu9   = 1026e-24;
sfpu9   = 741e-24;
scpu9   = sapu9-sfpu9;

% Pu-240
Npu0(1) = 0;
sfpu0   = 0.03e-24;
scpu0   = 289.5e-24;
sapu0   = scpu0+sfpu0;

% Pu-241
Npu1(1) = 0;
sfpu1   = 950e-24;
scpu1   = 362e-24;
sapu1   = sfpu1+scpu1;
htpu1   = 14*60*60*24*30*12;
lpu1    = log(2)/htpu1;

% Am-241
Nam1(1) = 0;
htam1   = 432.2*360*24*3600;
lam1    = log(2)/htam1;

% Pu-241
Npu2(1) = 0;
sfpu2   = 0.19e-24;
scpu2   = 30e-24;
sapu2   = sfpu2+scpu2;

% Total
tot(1)  = Nu8(1);
%}

%{
%======================== Looping ========================
i = 1;

while (t(i) ~= tf)
% U-238
    Nu8(i+1)    = (((n-1)*sau8*flux*dt)+1)*Nu8(i)/(1+sau8*flux*dt*n);
% U-239
    Nu9(i+1)    = (scu8*flux*dt*(n*Nu8(i+1)+(1-n)*Nu8(i))+(1-(1-n)*lu9*dt)*Nu9(i))/(1+lu9*dt*n);
% Np-239
    Nnp9(i+1)   = (lu9*dt*(n*Nu9(i+1)+(1-n)*Nu9(i))+(1-(1-n)*lnp9*dt)*Nnp9(i))/(1+lnp9*dt*n);
% Pu-239
    Npu9(i+1)   = (lnp9*dt*(n*Nnp9(i+1)+(1-n)*Nnp9(i))+(1-(1-n)*sapu9*flux*dt)*Npu9(i))/(1+sapu9*flux*dt*n);
% Pu-240
    Npu0(i+1)   = (scpu9*flux*dt*(n*Npu9(i+1)+(1-n)*Npu9(i))+(1-(1-n)*sapu0*flux*dt)*Npu0(i))/(1+sapu0*flux*dt*n);
% Pu-241
    Npu1(i+1)   = (scpu0*flux*dt*(n*Npu0(i+1)+(1-n)*Npu0(i))+(1-(1-n)*(lpu1+sapu1*flux)*dt)*Npu1(i))/(1+(lpu1+sapu1*flux)*dt*n);
% Am-241
    Nam1(i+1)   = (lpu1*dt*(n*Npu1(i+1)+(1-n)*Npu1(i))+(1-(1-n)*lam1*dt)*Nam1(i))/(1+lam1*dt*n); 
% Pu-242
    Npu2(i+1)   = (scpu1*flux*dt*(n*Npu1(i+1)+(1-n)*Npu1(i))+(1-(1-n)*sapu2*flux*dt)*Npu2(i))/(1+sapu2*flux*dt*n);
% Total
    tot(i+1) = Nu8(i+1)+Nu9(i+1)+Nnp9(i+1)+Npu9(i+1)+Npu0(i+1)+Npu1(i+1)+Nam1(i+1)+Npu2(i+1);
    i = i+1;
end
%}

%
%======================= Plotting Neutronik ===============================
figure(1)
plot(psirad); %Plot Neutronik Radial
grid minor
title('Grafik Neutronik Arah Radial');
xlabel('r (cell)');
ylabel('Fluks (cm^2/s)');

figure(2)
plot(psiax); %Plot Neutronik Aksial
grid minor
title('Grafik Neutronik Arah Aksial');
xlabel('h (cell)');
ylabel('Fluks (cm^2/s)');

figure(3)
mesh(powcellax); %Plot Daya per Satuan Luas
title('Grafik Daya per Satuan Luas');
xlabel('h (cell)');
ylabel('r (cell)');
zlabel('Daya/Luas (watt/m^2)');

%{
%======================= Plotting Burn-Up =================================
%detik*menit*jam*hari*bulan*tahun
t = t/(60*60*24*30*1*1);
figure (4);
plot (t,Nu8,'LineWidth',2)
grid minor
title('Jumlah Isotop U-238 terhadap Waktu');
xlabel('t (bulan)');
ylabel('Densitas Atom (atom/barn.cm)');
legend('U-238');

figure (5);
plot (t,Nu9,'LineWidth',2)
grid minor
title('Jumlah Isotop U-239 terhadap Waktu');
xlabel('t (bulan)');
ylabel('Densitas Atom (atom/barn.cm)');
legend('U-239');

figure (6);
plot (t,Nnp9,'LineWidth',2)
grid minor
title('Jumlah Isotop Np-239 terhadap Waktu');
xlabel('t (bulan)');
ylabel('Densitas Atom (atom/barn.cm)');
legend('Np-239');

figure (7);
plot (t,Npu9,'LineWidth',2)
grid minor
title('Jumlah Isotop Pu-239 terhadap Waktu');
xlabel('t (bulan)');
ylabel('Densitas Atom (atom/barn.cm)');
legend('Pu-239');

figure (8);
plot (t,Npu0,'LineWidth',2)
grid minor
title('Jumlah Isotop Pu-240 terhadap Waktu');
xlabel('t (bulan)');
ylabel('N (jumlah isotop)');
legend('Pu-240');

figure (9);
plot (t,Npu1,'LineWidth',2)
grid minor
title('Jumlah Isotop Pu-241 terhadap Waktu');
xlabel('t (bulan)');
ylabel('Densitas Atom (atom/barn.cm)');
legend('Pu-241');

figure (10);
plot (t,Nam1,'LineWidth',2)
grid minor
title('Jumlah Isotop Am-241 terhadap Waktu');
xlabel('t (bulan)');
ylabel('Densitas Atom (atom/barn.cm)');
legend('Am-241');

figure (11);
plot (t,Npu2,'LineWidth',2)
grid minor
title('Jumlah Isotop Pu-242 terhadap Waktu');
xlabel('t (bulan)');
ylabel('Densitas Atom (atom/barn.cm)');
legend('Pu-242');

figure (12);
plot (t,tot,'LineWidth',2)
grid minor
title('Jumlah Total Nuklida terhadap Waktu');
xlabel('t (bulan)');
ylabel('Densitas Atom (atom/barn.cm)');
legend('Total Nuklida');
%}
