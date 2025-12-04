%{
Script to generate scatter plots for Part 2 of BIOENG1005 Computer Project 2.

Data obtained from results of RF shimming.

Author: Paul Kullmann
%}


b1_color1 = [32 68 180]/180; % dark blue
b1_color2 = [85 128 203]/203;% lighter blue
sar_color1 = [220 97 39]/220;% orange
sar_color2 = [207 140 39]/207;% lighter orange

%% slices
slices = [45 50 55];
b1_cov = [0.2 0.05 0.11];
b1_ratio = [2.37 1.22 1.89]; 
sar_max = [17.7755 10.4341 7.9943];
sar_mean = [4.2781 3.2811 2.3785];

figure();
subplot(2,3,1)
grid on;
title("Individual Slices: B1+ Field")
yyaxis left;
scatter(slices, b1_cov, 'HandleVisibility', 'off', MarkerFaceColor=b1_color1); hold on;
plot(slices, b1_cov, 'DisplayName','B1+ COV', Color=[32 68 180 90]/180, LineWidth=3,LineStyle='-');
ylim([0 0.3]);
ylabel("B1+ COV");

yyaxis right;
scatter(slices, b1_ratio, 'HandleVisibility', 'off', MarkerFaceColor=sar_color1);
plot(slices, b1_ratio, 'DisplayName','B1+ Max/Min', Color=[220 97 39 110]/220, LineWidth=3,LineStyle='-');
ylim([0 10]);
ylabel("B1+ Max/Min Ratio");
xlabel("Slice Location");
xlim([42 58])
legend('Location','best');


subplot(2,3,4)
grid on;
title("Individual Slices: SAR");
yyaxis left;
scatter(slices, sar_max, 'HandleVisibility', 'off', MarkerFaceColor=b1_color1);
plot(slices, sar_max, 'DisplayName','Peak SAR', Color=[32 68 180 90]/180, LineWidth=3,LineStyle='-');
ylim([0 20]);
ylabel("Peak SAR (W/kg/10g)");

yyaxis right;
scatter(slices, sar_mean, 'HandleVisibility', 'off', MarkerFaceColor=sar_color1);
plot(slices, sar_mean, 'DisplayName','Mean SAR', Color=[220 97 39 110]/220, LineWidth=3,LineStyle='-');
ylabel("Mean SAR (W/kg/10g)");
xlim([42 58]);
xlabel("Slice Location");
legend;
ylim([0 6]);

%% 4-section slabs
slices = [mean([44 47]) mean([49 53]) mean([55 58])];
b1_cov = [0.21 0.13 0.23];
b1_ratio = [5.93 2.51 9.83]; 
sar_max = [10.1545 8.1765 5.8443];
sar_mean = [2.4315 2.7498 1.9393];

subplot(2,3,2)
grid on;
title("4-Slice Thick Sections: B1+ Field");
yyaxis left;
plot([44 47], [5.93 5.93], 'HandleVisibility','off', Color=[0 0 0 0.3], LineWidth=1.5);
scatter(slices, b1_cov, 'HandleVisibility', 'off', MarkerFaceColor=b1_color1); hold on;
plot(slices, b1_cov, 'DisplayName','B1+ COV', Color=[32 68 180 90]/180, LineWidth=3,LineStyle='-');
ylim([0 0.3]);
ylabel("B1+ COV");

yyaxis right;
scatter(slices, b1_ratio, 'HandleVisibility', 'off', MarkerFaceColor=sar_color1);
plot(slices, b1_ratio, 'DisplayName','B1+ Max/Min', Color=[220 97 39 110]/220, LineWidth=3,LineStyle='-');
ylim([0 10]);
ylabel("B1+ Max/Min Ratio");
xlabel("Center of Slab Location");
xlim([42 58])
legend('Location','best');


subplot(2,3,5)
grid on;
yyaxis left;
scatter(slices, sar_max, 'HandleVisibility', 'off', MarkerFaceColor=b1_color1);
plot(slices, sar_max, 'DisplayName','Peak SAR', Color=[32 68 180 90]/180, LineWidth=3,LineStyle='-');
ylim([0 20]);
ylabel("Peak SAR (W/kg/10g)");

yyaxis right;
scatter(slices, sar_mean, 'HandleVisibility', 'off', MarkerFaceColor=sar_color1);
plot(slices, sar_mean, 'DisplayName','Mean SAR', Color=[220 97 39 110]/220, LineWidth=3,LineStyle='-');
ylabel("Mean SAR (W/kg/10g)");
xlim([42 58]);
title("4-Slice Thick Sections: SAR");
xlabel("Center of Slab Location");
legend;
ylim([0 6]);



%% 8-section slabs
slices = [mean([42 49]) mean([47 54]) mean([51 58])];
b1_cov = [0.22 0.13 0.13];
b1_ratio = [4.45 2.48 3.35]; 
sar_max = [15.557 7.5612 8.8543];
sar_mean = [4.5059 2.4063 3.0187];

subplot(2,3,3)
grid on;
title("8-Slice Thick Sections: B1+ Field")
yyaxis left;
plot([42 49], [4.45 4.45], 'HandleVisibility','off', Color=[0 0 0 0.3], LineWidth=1.5);
scatter(slices, b1_cov, 'HandleVisibility', 'off', MarkerFaceColor=b1_color1); hold on;
plot(slices, b1_cov, 'DisplayName','B1+ COV', Color=[32 68 180 90]/180, LineWidth=3,LineStyle='-');
ylim([0 0.3]);
ylabel("B1+ COV");

yyaxis right;
scatter(slices, b1_ratio, 'HandleVisibility', 'off', MarkerFaceColor=sar_color1);
plot(slices, b1_ratio, 'DisplayName','B1+ Max/Min', Color=[220 97 39 110]/220, LineWidth=3,LineStyle='-');
ylim([0 10]);
ylabel("B1+ Max/Min Ratio");
xlim([42 58])
xlabel("Center of Slab Location");
legend('Location','best');


subplot(2,3,6);
grid on;
yyaxis left;
scatter(slices, sar_max, 'HandleVisibility', 'off', MarkerFaceColor=b1_color1);
plot(slices, sar_max, 'DisplayName','Peak SAR', Color=[32 68 180 90]/180, LineWidth=3,LineStyle='-');
ylim([0 20]);
ylabel("Peak SAR (W/kg/10g)");

yyaxis right;
scatter(slices, sar_mean, 'HandleVisibility', 'off', MarkerFaceColor=sar_color1);
plot(slices, sar_mean, 'DisplayName','Mean SAR', Color=[220 97 39 110]/220, LineWidth=3,LineStyle='-');
ylabel("Mean SAR (W/kg/10g)");


xlim([42 58]);
title("8-Slice Thick Sections: SAR");
xlabel("Center of Slab Location");
legend;
ylim([0 6]);
