
clear all; close all; clc;

rx1_ecef = [2176877.7800, 621906.9814, 5942863.6590];
rx2_ecef = [2176681.9919,622095.7950,5942910.4073];
rx3_ecef = [2176776.6563,622527.8168,5942839.3876];
rx4_ecef = [2176289.7709,622770.3445,5942992.6023];
rx5_ecef = [2176962.5171,622509.1930,5942777.2227];
C_ecef_to_enu = ecef_to_enu(69.2730616018062,15.9439788454853);
reference_enu = C_ecef_to_enu*rx1_ecef'
rx2_enu = C_ecef_to_enu*rx2_ecef';
rx3_enu = C_ecef_to_enu*rx3_ecef';
rx4_enu = C_ecef_to_enu*rx4_ecef';
rx5_enu = C_ecef_to_enu*rx5_ecef';

TX_ecef_check = lla2ecef([69.28007238, 16.00643461,381.98],'WGS84')
tx_enu = C_ecef_to_enu*TX_ecef_check';

pos2 = -(reference_enu-rx2_enu);
pos3 = -(reference_enu-rx3_enu);
pos4 = -(reference_enu-rx4_enu);
pos5 = -(reference_enu-rx5_enu);
postx = -(reference_enu-tx_enu);
diff2 = norm(pos2(1:2,1))
diff3 = norm(pos3(1:2,1))
diff4 = norm(pos4(1:2,1))
diff5 = norm(pos5(1:2,1))

diff23 = norm(rx2_enu(1:2,1)-rx3_enu(1:2,1));
diff45 = norm(rx4_enu(1:2,1)-rx4_enu(1:2,1));
diff24 = norm(rx2_enu(1:2,1)-rx4_enu(1:2,1));
diff34 = norm(rx3_enu(1:2,1)-rx4_enu(1:2,1));

% f = @(x,y) sqrt((x^2+y^2)) - sqrt((x-pos2(1))^2+(y-pos2(2))^2) - diff2;
% figure()
% fimplicit(f,[-200 1000 -500 1000])

f1 = @(x, y) abs(sqrt(x.^2 + y.^2) - sqrt((x - pos2(1)).^2 + (y - pos2(2)).^2)) - diff2;
f2 = @(x, y) abs(sqrt(x.^2 + y.^2) - sqrt((x - pos3(1)).^2 + (y - pos3(2)).^2)) - diff3;
f3 = @(x, y) abs(sqrt(x.^2 + y.^2) - sqrt((x - pos4(1)).^2 + (y - pos4(2)).^2)) - diff4;
f4 = @(x, y) abs(sqrt(x.^2 + y.^2) - sqrt((x - pos5(1)).^2 + (y - pos5(2)).^2)) - diff5;
f5 = @(x, y) abs(sqrt((x - pos2(1)).^2 + (y - pos2(2)).^2) - sqrt((x - pos3(1)).^2 + (y - pos3(2)).^2)) - diff23;
f6 = @(x, y) abs(sqrt((x - pos4(1)).^2 + (y - pos4(2)).^2) - sqrt((x - pos5(1)).^2 + (y - pos5(2)).^2)) - diff45;
f7 = @(x, y) abs(sqrt((x - pos2(1)).^2 + (y - pos2(2)).^2) - sqrt((x - pos4(1)).^2 + (y - pos4(2)).^2)) - diff23;
f8 = @(x, y) abs(sqrt((x - pos3(1)).^2 + (y - pos3(2)).^2) - sqrt((x - pos4(1)).^2 + (y - pos4(2)).^2)) - diff34;

[x, y] = meshgrid(-200:5:2700, -900:5:1400);  % Adjust the range as needed

% Evaluate the function on the meshgrid
z1 = f1(x, y);  % Apply the function element-wise over the grid
z2 = f2(x, y);
z3 = f3(x, y);
z4 = f4(x, y);
z5 = f5(x, y);
z6 = f6(x, y);
z7 = f7(x, y);
z8 = f8(x, y);

ztot = z1 + z2 + z3 + z4 + + z5 + z6 + z7 + z8;
% Create a 3D surface plot with color representing the function value
% Create a 2D filled contour plot
figure()
contourf(x, y, ztot, 20);  % 20 levels of contour
hold on
colorbar;  % Add a color bar to show the scale of values
plot(0,0,"r*",'MarkerSize', 12)
plot(pos2(1),pos2(2),"r*",'MarkerSize', 12);
plot(pos3(1),pos3(2),"r*",'MarkerSize', 12);
plot(pos4(1),pos4(2),"r*",'MarkerSize', 12);
plot(pos5(1),pos5(2),"r*",'MarkerSize', 12);
plot(postx(1),postx(2),"b*",'MarkerSize', 12);
title('Contour Plot of Potential TX Locations in ENU Relative to RX1');
axis equal
xlabel('X [m]');
ylabel('Y [m]');
hold off

figure()
plot(pos2(1),pos2(2),"r*",'MarkerSize', 12);
hold on
plot(pos3(1),pos3(2),"r*",'MarkerSize', 12);
plot(pos4(1),pos4(2),"r*",'MarkerSize', 12);
plot(pos5(1),pos5(2),"r*",'MarkerSize', 12);
plot(postx(1),postx(2),"b*",'MarkerSize', 12);
plot(0,0,"r*",'MarkerSize', 12)
xlim([-200,2700])
ylim([-500,1000])
title('Locations of RXs and TX in ENU Relative to RX1');
xlabel('X [m]');
ylabel('Y [m]');
grid minor
axis equal
hold off

% 
% [min_value, linear_index] = max(ztot(:));
% [row, col] = ind2sub(size(ztot), linear_index);
