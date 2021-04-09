%% входные данные
% объект в сферических координатах
psi = 0.5 * pi / 180;  lyam = 1 * pi / 180;
q0 = [psi, lyam];
% координаты спутников и расстояния до них
CO1 = [8757.75; -2183.55;  4305.11];
CO2 = [9061.25;  2598.27; -3338.07];
dist = zeros(2, 2);
dist(:, 1) = [5506.08; 4960];
dist(:, 2) = [5594.37; 4889.03];

%% алгоритм
% объект в прямоугольных координатах
r_earth = 6300;
x = @(psi, lyam)r_earth * cos(psi) * cos(lyam);
y = @(psi, lyam)r_earth * cos(psi) * sin(lyam);
z = @(psi, lyam)r_earth * sin(psi);
 
% счислимое расстояние до n-го спутника
dist_sp = @(psi, lyam, sp)sqrt((x(psi, lyam) - sp(1)).^2 + ...
                               (y(psi, lyam) - sp(2)).^2 + ...
                               (z(psi, lyam) - sp(3)).^2);

% част. производные
ddist_psi  = @(psi, lyam, sp)r_earth ./ dist_sp(psi, lyam, sp) .* (sp(1) * cos(lyam) * sin(psi) ...
                - sp(2) * sin(lyam) * sin(psi) - sp(3) * cos(psi));
ddist_lyam = @(psi, lyam, sp)r_earth ./ dist_sp(psi, lyam, sp) .* cos(psi) .* ...
                 (sp(1) * sin(lyam) - sp(2) * cos(lyam));
 
% матрица В
B = @(psi, lyam, sp1, sp2)[ddist_psi(psi, lyam, sp1) ddist_lyam(psi, lyam, sp1); ...
                           ddist_psi(psi, lyam, sp2) ddist_lyam(psi, lyam, sp2)]; 

% уточнение расстояния до спутников
delta_d = @(psi, lyam, sp1, sp2, dist)[dist(1) - dist_sp(psi, lyam, sp1); ...
                                       dist(2) - dist_sp(psi, lyam, sp2)];
% уточнение координат
delta_q = @(psi, lyam, sp1, sp2, dist)B(psi, lyam, sp1, sp2) \ delta_d(psi, lyam, sp1, sp2, dist);
 
% новые координаты
new_q = @(psi, lyam, sp1, sp2, dist)[psi; lyam] + delta_q(psi, lyam, sp1, sp2, dist);

%% результаты первого приближения
for i = 1:3
    q1 = new_q(psi, lyam, CO1, CO2, dist(:, 1));
    psi = q1(1);
    lyam = q1(2);
end

%% результаты второго приближения
for i = 1:3
    q2 = new_q(psi, lyam, CO1, CO2, dist(:, 2));
    psi = q2(1);
    lyam = q2(2);
end

%% задача идентификации
disp(['Изначальные координаты:    [', num2str(q0(1) * 180 / pi, '%10.3f'), ', ', num2str(q0(2) * 180 / pi, '%10.3f'), ']']);
disp(['Координаты в момент t1:    [', num2str(q1(1) * 180 / pi, '%10.3f'), ', ', num2str(q1(2) * 180 / pi, '%10.3f'), ']']);
disp(['Координаты в момент t2:    [', num2str(q2(1) * 180 / pi, '%10.3f'), ', ', num2str(q2(2) * 180 / pi, '%10.3f'), ']']);
% матрица перехода
S0 = [q0', q1]; S1 = [q1, q2];
F = S1 / S0;
disp('Матрица перехода: ');
disp(F);