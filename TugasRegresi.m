%Contoh Soal
% Data tabel
x = [0.0, 1.0, 2.0, 2.5, 3.0];
y = [2.9, 3.7, 4.1, 4.4, 5.0];

% Mencari nilai a dan b
n = length(x);
xy = x .* y;
x2 = x .^ 2;

b = (n * sum(xy) - sum(x) * sum(y)) / (n * sum(x2) - sum(x) ^ 2);
a = mean(y) - b * mean(x);

% Persamaan garis regresi
regression_line = a + b * x;

fprintf('\nContoh Soal\n')
% Menampilkan hasil
fprintf('Persamaan garis regresi: y = %.2f + %.2fx\n', a, b);

% Menampilkan tabel
fprintf('No\t\tx\t\ty\t\txy\t\tx^2\n');
for i = 1:n
    fprintf('%d\t %.1f\t %.1f\t %.1f\t %.1f\n', i, x(i), y(i), xy(i), x2(i));
end

% Menampilkan hasil tambahan
jumx = sum(x);
jumxy = sum(x .* y);
fprintf('Jumlah kolom x: %.4f\n', jumx);
fprintf('Jumlah kolom xy: %.4f\n', jumxy);

% Menampilkan grafik hasil regresi
scatter(x, y, 'Marker', 'o', 'DisplayName', 'Data Asli');
hold on;
plot(x, regression_line, 'r', 'DisplayName', 'Regresi Linear');
xlabel('x');
ylabel('y');
legend('Location', 'northwest');
title('Regresi Linear');
grid on;
hold off;



%No.1
% Data tabel
x = [0.5, 1.0, 1.5, 2.0, 2.5];
y = [0.49, 1.60, 3.36, 6.44, 10.16];

% Mencari nilai a dan b
lnx = log(x);
lny = log(y);

n = length(x);
lnxy = lnx .* lny;
lnx2 = lnx .^ 2;

b = (n * sum(lnxy) - sum(lnx) * sum(lny)) / (n * sum(lnx2) - sum(lnx) ^ 2);
a = exp(mean(lny) - b * mean(lnx));

fprintf('\n\nNomer 1\n')

% Persamaan garis regresi
fprintf('Persamaan garis regresi: y = %.2fx^%.2f\n', a, b);

% Menampilkan tabel
fprintf('No\t\tx\t\ty\t\txy\t\tx^2\n');
for i = 1:n
    fprintf('%d\t %.2f\t %.2f\t %.2f\t %.2f\n', i, x(i), y(i), x(i)*y(i), x(i)^2);
end

% Menampilkan hasil tambahan
jumx = sum(x);
jumxy = sum(x .* y);
fprintf('Jumlah kolom x: %.4f\n', jumx);
fprintf('Jumlah kolom xy: %.4f\n', jumxy);

% Menampilkan grafik hasil regresi
x_regresi = 0:0.1:3; % Rentang x untuk plot hasil regresi
y_regresi = a * exp(b * log(x_regresi));

figure;
plot(x, y, 'o', 'DisplayName', 'Data Asli');
hold on;
plot(x_regresi, y_regresi, 'r', 'DisplayName', 'Regresi Linear');
xlabel('x');
ylabel('y');
legend('Location', 'northwest');
title('Regresi Linear');
grid on;
hold off;


%No. 2
% Data yang diberikan
t = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5];
gamma = [1.000, 0.994, 0.990, 0.985, 0.979, 0.977, 0.972, 0.969, 0.967, 0.960, 0.956, 0.952];

% Nilai a dan b yang diberikan
a = 0.9984;
b = 0.0086;

fprintf('\n\nNomer 2\n')

% Menampilkan hasil
fprintf('Nilai a = %.4f\n', a);
fprintf('Nilai b = %.4f\n', b);

% Mencetak tabel tambahan
fprintf('No\t\tt\t\tγ\t\ttγ\t\tt^2\n');
for i = 1:length(t)
    fprintf('%d\t%.2f\t%.3f\t%.3f\t%.2f\n', i, t(i), gamma(i), t(i) * gamma(i), t(i)^2);
end

% Menghitung jumlah kolom t dan jumlah kolom tγ
sum_t = sum(t);
sum_t_gamma = sum(t .* gamma);

fprintf('Jumlah kolom t: %.4f\n', sum_t);
fprintf('Jumlah kolom tγ: %.4f\n', sum_t_gamma);

% Menampilkan perubahan pada teks Persamaan garis regresi
fprintf('Persamaan garis regresi: γ(t) = %.4f * exp^(-%.4f * t)\n', a, b);

% Mencetak grafik hasil regresi
t_fit = linspace(0, 6, 100); % Rentang nilai t untuk plot hasil regresi
gamma_fit = a * exp(-b * t_fit);

figure;
plot(t, gamma, 'o', t_fit, gamma_fit);
xlabel('t');
ylabel('\gamma(t)');
legend('Data', 'Regresi');

