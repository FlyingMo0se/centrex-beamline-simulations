%% Import data for plotting %%
filename = 'Data/loop_over_voltage_06_27_2018_J=2_d=1.75in_v=200_StoL=95cm_N=1e7_RK4.txt';
delimiter = '\t';
header_lines = 1;
data = importdata(filename, delimiter, header_lines);

%% Plot gain vs lens length
figure
hold on
errorbar(data.data(:,5), data.data(:,6),data.data(:,7),'x')
title('Efficiency vs lens voltage,S to L = 0.95 m, v = 200 m/s, d = 1.75 in')
xlabel(data.colheaders(5))
ylabel(data.colheaders(6))
%saveas(gcf,'./Plots/Efficiency vs lens length 6_21_2018, d = 1.75 in.pdf')

% %% Plot fraction hitting field plates
% figure
% errorbar(data.data(:,1), data.data(:,8),data.data(:,9),'-r')
% title('f_{FP} vs lens length')
% xlabel(data.colheaders(1))
% ylabel(data.colheaders(8))
% saveas(gcf,'f_fp_vs_length_zofzd=0.75cm_cell_to_zoz=1cm_N=1e7_2.pdf')

%% Plot gain vs electrode length
figure
hold on

errorbar(data.data(:,5), data.data(:,8),data.data(:,9),'x')
title('Gain vs source to lens distance, v_z = 200 m/s, d = 1.75 in')
xlabel(data.colheaders(5))
ylabel(data.colheaders(8))
%saveas(gcf,'./Plots/Gain vs lens length 6_21_2018, d = 1.75 in.pdf')
