% Name = Naveen Srinivasa
clear;
close all;
clc;

seed = 16+64+64+256+256+16+16;   % a = 16, e = 64, i = 256 seed = 688
seed
rng(seed,'twister');
% Given Values
m =1000;  % Number of trials
C=10;  % True Channel Coefficient
sigma_n = 20;  % Standard deviation of the noise
confidence = 0.9; % 90% confidence
alpha = 1-confidence;
% Performing the simulation separately for n = 10 and n = 100

% n = 10
n_10 = 10;

count_case1_10 = 0;  % Number of cases where C lies in the confidence interval for case 1
count_case2_10 = 0;  % Number of cases where C lies in the confidence interval for case 2
count_case3_10 = 0;  % Number of cases where C lies in the confidence interval for case 3
sample_mean_list_10 = [0,0,0,0,0,0,0,0,0,0];  % Sample means from first 10 trials
lower_conf_1_10 = [0,0,0,0,0,0,0,0,0,0];  % Lower limit of the confidence interval for first 10 trials for case 1
upper_conf_1_10 = [0,0,0,0,0,0,0,0,0,0];  % Upper limit of the confidence interval for first 10 trials for case 1
lower_conf_2_10 = [0,0,0,0,0,0,0,0,0,0];  % Lower limit of the confidence interval for first 10 trials for case 2
upper_conf_2_10 = [0,0,0,0,0,0,0,0,0,0];  % Upper limit of the confidence interval for first 10 trials for case 2
lower_conf_3_10 = [0,0,0,0,0,0,0,0,0,0];  % Lower limit of the confidence interval for first 10 trials for case 3
upper_conf_3_10 = [0,0,0,0,0,0,0,0,0,0];  % Upper limit of the confidence interval for first 10 trials for case 3

%Loop for 1000 trials
for i = 1:m
    N = sigma_n*randn(1,n_10);   % Noise
    X = C+N;  % Real valued Received Sample
    sample_mean = sum(X)/n_10;  % Sample Mean for the generated Data points

    %case 1 : Estimating the confidence interval using known variance of X
    sigma_1 = sigma_n; % For case 1, The variance of the channel is the same as that of the noise
    delta_1_10 = norminv((1-alpha/2),0,sigma_1/sqrt(n_10)); % Error
    y_value_1_10 = delta_1_10 * sqrt(n_10)/sigma_1;  
    if ((sample_mean-delta_1_10<=C) && (C <= sample_mean+delta_1_10))
        count_case1_10 = count_case1_10+1;  %Counting the number of times when C lies in the confidence interval
    end
    

    %case 2 : Estimating the confidence interval using Sample variance without assuming that the measurements are gaussian
    sample_variance = (sum((X-sample_mean).^2))/(n_10-1);
    sigma_2 = sqrt(sample_variance);
    y_value_2_10 = y_value_1_10;
    delta_2_10 = y_value_2_10*sigma_2/sqrt(n_10);  % Error
    if ((sample_mean-delta_2_10<=C) && (C <= sample_mean+delta_2_10))
        count_case2_10 = count_case2_10+1;  %Counting the number of times when C lies in the confidence interval
    end
   

    %case 3 : Estimating the confidence interval using Sample variance assuming that the measurements are gaussian
    y_value_3_10 = tinv(1-alpha/2,n_10-1);
    sigma_3 = sigma_2;
    delta_3_10 = y_value_3_10*sigma_3/sqrt(n_10);  
    if ((sample_mean-delta_3_10<=C) && (C <= sample_mean+delta_3_10))
        count_case3_10 = count_case3_10+1;  %Counting the number of times when C lies in the confidence interval
    end
    
    if (i<=10)
        % Storing the 1st 10 values of the sample mean and the confidence
        % interval bounds generated for the 3 cases
        sample_mean_list_10(i) = sample_mean;
        lower_conf_1_10(i) = sample_mean - delta_1_10;
        upper_conf_1_10(i) = sample_mean + delta_1_10;
        lower_conf_2_10(i) = sample_mean - delta_2_10;
        upper_conf_2_10(i) = sample_mean + delta_2_10;
        lower_conf_3_10(i) = sample_mean - delta_3_10;
        upper_conf_3_10(i) = sample_mean + delta_3_10;
    end

    
end


figure;
x_data = linspace(1,10,10);
x_line = linspace(0,12,12);
y_line = ones(1, 12)*C;
plot(x_line, y_line, 'r'); hold on;
err1 = upper_conf_1_10 - sample_mean_list_10;
errorbar(x_data, sample_mean_list_10, err1, 'k.'); % Plots error bars at the specific data points given by x_data and y_data.
plot(x_data, sample_mean_list_10, 'b.', 'MarkerSize', 15); hold on;% Add data points to the plot.


xlim([0 11]);
ylim([-25 35]);
title('Case 1 (n=10)');
xlabel('Number of Trials');



figure;
x_line = linspace(0,12,12);
y_line = ones(1, 12)*C;
plot(x_line, y_line, 'r'); hold on;
err2 = upper_conf_2_10 - sample_mean_list_10;
errorbar(x_data, sample_mean_list_10, err2, 'k.'); % Plots error bars at the specific data points given by x_data and y_data.
plot(x_data, sample_mean_list_10, 'b.', 'MarkerSize', 15); hold on;% Add data points to the plot.


xlim([0 11]);
ylim([-25 35]);
title('Case 2 (n=10)');
xlabel('Number of Trials');

legend('Channel Coefficient', 'Confidence Interval', 'Sample Mean');
hold off

figure;
x_line = linspace(0,12,12);
y_line = ones(1, 12)*C;
plot(x_line, y_line, 'r'); hold on;
err3 = upper_conf_3_10 - sample_mean_list_10;
errorbar(x_data, sample_mean_list_10, err3, 'k.'); % Plots error bars at the specific data points given by x_data and y_data.
plot(x_data, sample_mean_list_10, 'b.', 'MarkerSize', 15); hold on;% Add data points to the plot.


xlim([0 11]);
ylim([-25 35]);
title('Case 3 (n=10)');
xlabel('Number of Trials');

legend('Channel Coefficient', 'Confidence Interval', 'Sample Mean');
hold off

% n = 100
n_100 = 100;

count_case1_100 = 0;  % Number of cases where C lies in the confidence interval for case 1
count_case2_100 = 0;  % Number of cases where C lies in the confidence interval for case 2
count_case3_100 = 0;  % Number of cases where C lies in the confidence interval for case 3
sample_mean_list_100 = [0,0,0,0,0,0,0,0,0,0];  % Sample means from first 10 trials
lower_conf_1_100 = [0,0,0,0,0,0,0,0,0,0];  % Lower limit of the confidence interval for first 10 trials for case 1
upper_conf_1_100 = [0,0,0,0,0,0,0,0,0,0];  % Upper limit of the confidence interval for first 10 trials for case 1
lower_conf_2_100 = [0,0,0,0,0,0,0,0,0,0];  % Lower limit of the confidence interval for first 10 trials for case 2
upper_conf_2_100 = [0,0,0,0,0,0,0,0,0,0];  % Upper limit of the confidence interval for first 10 trials for case 2
lower_conf_3_100 = [0,0,0,0,0,0,0,0,0,0];  % Lower limit of the confidence interval for first 10 trials for case 3
upper_conf_3_100 = [0,0,0,0,0,0,0,0,0,0];  % Upper limit of the confidence interval for first 10 trials for case 3
%Loop for 1000 trials
for i = 1:m
    N = sigma_n*randn(1,n_100);   % Noise
    X = C+N;  % Real valued Received Sample
    sample_mean = sum(X)/n_100;  % Sample Mean for the generated Data points

    %case 1 : Estimating the confidence interval using known variance of X
    sigma_1 = sigma_n; % For case 1, The variance of the channel is the same as that of the noise
    delta_1_100 = norminv((1-alpha/2),0,sigma_1/sqrt(n_100)); % Error
    y_value_1_100 = delta_1_100 * sqrt(n_100)/sigma_1;  
    if ((sample_mean-delta_1_100<=C) && (C <= sample_mean+delta_1_100))
        count_case1_100 = count_case1_100+1;  %Counting the number of times when C lies in the confidence interval
    end
    

    %case 2 : Estimating the confidence interval using Sample variance without assuming that the measurements are gaussian
    sample_variance = (sum((X-sample_mean).^2))/(n_100-1);
    sigma_2 = sqrt(sample_variance);
    y_value_2_100 = y_value_1_100;
    delta_2_100 = y_value_2_100*sigma_2/sqrt(n_100);  % Error
    if ((sample_mean-delta_2_100<=C) && (C <= sample_mean+delta_2_100))
        count_case2_100 = count_case2_100+1;  %Counting the number of times when C lies in the confidence interval
    end
   

    %case 3 : Estimating the confidence interval using Sample variance assuming that the measurements are gaussian
    y_value_3_100 = tinv(1-alpha/2,n_100-1);
    sigma_3 = sigma_2;
    delta_3_100 = y_value_3_100*sigma_3/sqrt(n_100);  

    if ((sample_mean-delta_3_100<=C) && (C <= sample_mean+delta_3_100))
        count_case3_100 = count_case3_100+1;  %Counting the number of times when C lies in the confidence interval
    end
    
    if (i<=10)
        % Storing the 1st 10 values of the sample mean and the confidence
        % interval bounds generated for the 3 cases
        sample_mean_list_100(i) = sample_mean;
        lower_conf_1_100(i) = sample_mean - delta_1_10;
        upper_conf_1_100(i) = sample_mean + delta_1_10;
        lower_conf_2_100(i) = sample_mean - delta_2_10;
        upper_conf_2_100(i) = sample_mean + delta_2_10;
        lower_conf_3_100(i) = sample_mean - delta_3_10;
        upper_conf_3_100(i) = sample_mean + delta_3_10;
    end    
end