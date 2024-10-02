clc
clear all
format longG

% =============================
% Langkah 1: Input Data Pengamatan
% =============================

% Data observasi ketinggian air laut
data = readtable('../Data Pasang Surut BIG/datapasut18062015_16072015.txt', 'Delimiter', '\t', 'HeaderLines', 1);  % Delimiter diatur sebagai spasi, HeaderLines untuk melewatkan header
    
t = 0:1:height(data)-1; % Waktu dalam jam (25 data poin, dari jam 0 sampai 24)
t=transpose(t);

h_obs = table2array(data);  % Data observasi ketinggian air laut

% =============================
% Langkah 2: Definisi Komponen Pasut
% =============================

% Periode untuk 9 komponen pasut (dalam jam, dikonversi ke detik)
omega = [0.5059; 0.5236; 0.4964; 0.5250; 0.2625; 0.2434; 0.2611; 1.0117; 1.0295];
konstanta_nama = {'M2'; 'S2'; 'N2'; 'K2'; 'K1'; 'O1'; 'P1'; 'M4'; 'MS4'}; % Nama komponen

% =============================
% Langkah 3: Buat Matriks Desain
% =============================

% Jumlah data
n = length(t);
u = 2*height(omega)+1;

% Tambahkan kolom untuk cos(w_i * t) dan sin(w_i * t) untuk setiap komponen
for i = 1:n %t sebanyak 15 piantan (360 data)
    A(i,1) = 1; %kolom 1 pada matriks A
    for j = 1:length(konstanta_nama) %sebanyak konstanta
        A(i,2*j) = cos(omega(j,1)*i); %cos dimulai dr kolom ke 2,4,6,...,18
        A(i,2*j+1) = -sin(omega(j,1)*i); %sin dimulai dr kolom 3,5,7,...,19
    end
end

% =============================
% Langkah 4: Menghitung Parameter Least Square
% =============================
X = inv(A'*A)*A'*h_obs;
V = A*X-h_obs;
variansi_aposteori = (V'*V)/(n-u);
standar_deviasi = sqrt(variansi_aposteori);

% =============================
% Langkah 5: Menghitung Konstanta Pasut
% =============================
% Ekstrak H0
H0_est = X(1);

% Ekstrak C_i dan D_i
C = X(2:2:end);
D = X(3:2:end);

% Hitung Amplitudo dan Fase
A_est = sqrt(C.^2 + D.^2);
phi_est = atan2(D, C); % fase dalam radian

% Sesuaikan fase agar berada di antara 0 dan 2pi
phi_est(phi_est < 0) = phi_est(phi_est < 0) + 2*pi;

% =============================
% Langkah 6: Tampilkan Hasil
% =============================

fprintf('Estimasi Mean Sea Level (H0): %.4f meter\n\n', H0_est);

fprintf('Estimasi Komponen Pasang Surut:\n');
fprintf('%-5s %-12s %-12s %-12s %-12s %-10s\n', 'No.', 'Komponen','Amplitudo(m)', 'Fase(rad)', 'Fase(deg)');
for i = 1:length(konstanta_nama)
    fprintf('%-5d %-12s %-12.4f %-12.4f %-10.2f\n', i, konstanta_nama{i}, A_est(i), phi_est(i), rad2deg(phi_est(i)));
end

% =============================
% Langkah 6: Menentukan Jenis Pasut
% =============================
Formzahl = (A_est(5)+A_est(6))/(A_est(1)+A_est(2))

if 0<Formzahl<0.25
    fprintf("nilai formzahl sebesar %.7f (Semi Diurnal)\n\n",Formzahl)
    else if 0.25<Formzahl<1.5
    fprintf("nilai formzahl sebesar %.7f (Campuran Berganda)\n\n",Formzahl)
    else if 1.5<Formzahl<3
    fprintf("nilai formzahl sebesar %.7f (Campuran Tunggal)\n\n",Formzahl)
    else Formzahl>3
    fprintf("nilai formzahl sebesar %.7f (Diurnal Murni)\n\n",Formzahl)
    end 
    end
end


% =============================
% Langkah 7: Rekonstruksi Model
% =============================
% Rekonstruksi model ketinggian air laut menggunakan parameter yang diestimasi
h_model_est = H0_est * ones(size(t));
for i = 1:length(konstanta_nama)
    h_model_est = h_model_est + A_est(i) * cos(omega(i) * t + phi_est(i));
end

% =============================
% Langkah 8: Perhitungan Kedudukan Muka Air
% =============================

HHWL = H0_est+(A_est(1)+A_est(2)+A_est(4)+A_est(5)+A_est(6)+A_est(7));
fprintf("didapatkan nilai HHWL (High Highest Water Level) : %.5f meter\n", HHWL)

MHWL = H0_est+(A_est(1)+A_est(5)+A_est(6));
fprintf("didapatkan nilai MHWL (Mean Highest Water Level) : %.5f meter\n", MHWL)

MSL = H0_est;
fprintf("didapatkan nilai MSL (Mean Sea Level) : %.5f meter\n", MSL)

MLWL = H0_est-(A_est(1)+A_est(5)+A_est(6));
fprintf("didapatkan nilai MLWL (Mean Lowest Water Level) : %.5f meter\n", MLWL)

CDL = H0_est-(A_est(1)+A_est(2)+A_est(5)+A_est(6));
fprintf("didapatkan nilai CDL (Chart Datum Level) : %.5f meter\n", CDL)

LLWL = H0_est-(A_est(1)+A_est(2)+A_est(4)+A_est(5)+A_est(6)+A_est(7));
fprintf("didapatkan nilai LLWL (Low Lowest Water Level) : %.5f meter\n", LLWL)

LAT = H0_est-(A_est(1)+A_est(2)+A_est(3)+A_est(4)+A_est(5)+A_est(6)+A_est(7)+A_est(8)+A_est(9));
fprintf("didapatkan nilai LAT (Lowest Astronomical Tide) : %.5f meter\n\n\n", LAT)

% =============================
% Langkah 9: Visualisasi Hasil
% =============================

% Plot data observasi vs model
figure;
plot(t, h_obs, '-', 'DisplayName', 'Data Observasi');
hold on;
plot(t, V, '-', 'DisplayName', 'Koreksi');
hold on;
plot(t, h_model_est, '-', 'DisplayName', 'Model Pasang Surut');

% Menambahkan garis-garis horizontal untuk nilai-nilai kedudukan muka air
yline(HHWL, '--r', 'DisplayName', 'HHWL');
yline(MHWL, '--g', 'DisplayName', 'MHWL');
yline(MSL, '--b', 'DisplayName', 'MSL');
yline(MLWL, '--c', 'DisplayName', 'MLWL');
yline(CDL, '--y', 'DisplayName', 'CDL');
yline(LLWL, '--m', 'DisplayName', 'LLWL');
yline(LAT, '--k', 'DisplayName', 'LAT');

% Tambahkan label sumbu dan judul
xlabel('Waktu (jam)');
ylabel('Ketinggian Air Laut (m)');
title('Perbandingan Ketinggian Air Laut Observasi dan Model');

% Dapatkan batas sumbu x
xLimits = get(gca, 'XLim');
xticks(0:24:xLimits(2));

% Tampilkan legenda
legend("show", 'Orientation', 'horizontal')
grid on;
hold off;

RMSE = rmse(h_model_est,h_obs)