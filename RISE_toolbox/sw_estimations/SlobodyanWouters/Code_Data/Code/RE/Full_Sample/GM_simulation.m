function GM_simulation(TT,RR,Q,mus,lgy,order)

% diary Great_Moderation_T13d_5Mod.txt
% diary on
% pack;

RUNS = 5000;
num_shock = size(Q,1);

smpl = size(TT,3);

SState = (eye(size(TT,1)) - TT(:,:,1)) \ mus(:,1);

y = zeros(size(TT,1),smpl+1,RUNS);
y(:,1,:) = kron(SState,ones(1,RUNS));

randn('state',sum(100 * clock));
ex = randn(num_shock,RUNS,smpl);

for t = 1:smpl
    y(:,t+1,:) = kron(mus(:,t),ones(1,RUNS)) + TT(:,:,t) * squeeze(y(:,t,:)) + ...
        RR(:,:,t) * chol(Q) * squeeze(ex(:,:,t));
end

tmp = deblank(lgy(order,:));

pinf_pos = strmatch('pinfobs',tmp,'exact');
dy_pos   = strmatch('dy',tmp,'exact');
dc_pos   = strmatch('dc',tmp,'exact');
dinve_pos   = strmatch('dinve',tmp,'exact');
dw_pos   = strmatch('dw',tmp,'exact');
y_pos    = strmatch('y',tmp,'exact');
c_pos    = strmatch('c',tmp,'exact');

start_smpl = 71;
start_GM = 77;

% looking at real data's volatilities
load usmodel_data_SW_data_april2009 dy dc pinfobs;

yobs = cumsum(dy(start_smpl:end));
cobs = cumsum(dc(start_smpl:end));
yobs_f = HP(yobs,1600,0);
cobs_f = HP(cobs,1600,0);
data_std_pinf_bGM = std(pinfobs(start_smpl:(start_smpl+start_GM-1)));
data_std_pinf_aGM = std(pinfobs(start_smpl+start_GM:end));
data_std_dy_bGM = std(dy(start_smpl:(start_smpl+start_GM-1)));
data_std_dy_aGM = std(dy(start_smpl+start_GM:end));
data_std_dc_bGM = std(dc(start_smpl:(start_smpl+start_GM-1)));
data_std_dc_aGM = std(dc(start_smpl+start_GM:end));
data_std_yobs_f_bGM = std(yobs_f(1:start_GM-1))
data_std_yobs_f_aGM = std(yobs_f(start_GM:end))
data_std_cobs_f_bGM = std(cobs_f(1:start_GM-1))
data_std_cobs_f_aGM = std(cobs_f(start_GM:end))

data_pinf_GM = (data_std_pinf_bGM - data_std_pinf_aGM) / data_std_pinf_bGM;
data_dy_GM = (data_std_dy_bGM - data_std_dy_aGM) / data_std_dy_bGM;
data_dc_GM = (data_std_dc_bGM - data_std_dc_aGM) / data_std_dc_bGM;
data_yobs_f_GM = (data_std_yobs_f_bGM - data_std_yobs_f_aGM) / data_std_yobs_f_bGM;
data_cobs_f_GM = (data_std_cobs_f_bGM - data_std_cobs_f_aGM) / data_std_cobs_f_bGM;

std_y_aGM = zeros(RUNS,1);
std_y_bGM = zeros(RUNS,1);
std_c_aGM = zeros(RUNS,1);
std_c_bGM = zeros(RUNS,1);

for i = 1:RUNS
%     y_f = BK(squeeze(y(y_pos,:,i))',6,32,12);
    y_f = HP(squeeze(y(y_pos,:,i))',1600,0);
    std_y_bGM(i) = std(y_f(1:start_GM-1));
    std_y_aGM(i) = std(y_f(start_GM:end));    
%     c_f = BK(squeeze(y(c_pos,:,i))',6,32,12);
    c_f = HP(squeeze(y(c_pos,:,i))',1600,0);
    std_c_bGM(i) = std(c_f(1:start_GM-1));
    std_c_aGM(i) = std(c_f(start_GM:end));    
end

% variables vs 1984:Q1 to 2004:Q4. First observation is excluded - it's REE
% value; other four are the KF training sample
% Over-ruled, start from the beginning
mean_pinf_bGM = squeeze(mean(y(pinf_pos,1:start_GM-1,:)));
mean_pinf_aGM = squeeze(mean(y(pinf_pos,start_GM:end,:)));
mean_dy_bGM = squeeze(mean(y(dy_pos,1:start_GM-1,:)));
mean_dy_aGM = squeeze(mean(y(dy_pos,start_GM:end,:)));
mean_dc_bGM = squeeze(mean(y(dc_pos,1:start_GM-1,:)));
mean_dc_aGM = squeeze(mean(y(dc_pos,start_GM:end,:)));
mean_dinve_bGM = squeeze(mean(y(dinve_pos,1:start_GM-1,:)));
mean_dinve_aGM = squeeze(mean(y(dinve_pos,start_GM:end,:)));
mean_dw_bGM = squeeze(mean(y(dw_pos,1:start_GM-1,:)));
mean_dw_aGM = squeeze(mean(y(dw_pos,start_GM:end,:)));

std_pinf_bGM = squeeze(std(y(pinf_pos,1:start_GM-1,:)));
std_pinf_aGM = squeeze(std(y(pinf_pos,start_GM:end,:)));
std_dy_bGM = squeeze(std(y(dy_pos,1:start_GM-1,:)));
std_dy_aGM = squeeze(std(y(dy_pos,start_GM:end,:)));
std_dc_bGM = squeeze(std(y(dc_pos,1:start_GM-1,:)));
std_dc_aGM = squeeze(std(y(dc_pos,start_GM:end,:)));
std_dinve_bGM = squeeze(std(y(dinve_pos,1:start_GM-1,:)));
std_dinve_aGM = squeeze(std(y(dinve_pos,start_GM:end,:)));
std_dw_bGM = squeeze(std(y(dw_pos,1:start_GM-1,:)));
std_dw_aGM = squeeze(std(y(dw_pos,start_GM:end,:)));

% Percentages of GM explained
GM_pinf_ = (std_pinf_bGM - std_pinf_aGM) ./ std_pinf_bGM / data_pinf_GM * 100;
GM_pinf = mean(GM_pinf_);
GM_dy_ = (std_dy_bGM - std_dy_aGM) ./ std_dy_bGM / data_dy_GM * 100;
GM_dy = mean(GM_dy_);
GM_dc_ = (std_dc_bGM - std_dc_aGM) ./ std_dc_bGM / data_dc_GM * 100;
GM_dc = mean(GM_dc_);
GM_yobs_ = (std_y_bGM - std_y_aGM) ./ std_y_bGM / data_yobs_f_GM * 100;
GM_yobs = mean(GM_yobs_);
GM_cobs_ = (std_c_bGM - std_c_aGM) ./ std_c_bGM / data_cobs_f_GM * 100;
GM_cobs = mean(GM_cobs_);

% Percentage of GMs exceeding the data
count_GM_pinf = length(find(GM_pinf_ > 100)) / RUNS * 100;
count_GM_dy   = length(find(GM_dy_ > 100))   / RUNS * 100;
count_GM_dc   = length(find(GM_dc_ > 100))   / RUNS * 100;
count_GM_yobs = length(find(GM_yobs_ > 100)) / RUNS * 100;
count_GM_cobs = length(find(GM_cobs_ > 100)) / RUNS * 100;

disp('RESULTS OF THE SIMULATION');
fprintf('  Pre-84:Mean Post-84:Mean   Pre-84:Std Post-84:Std \n');
fprintf('pinf  %10.2f %10.2f %10.2f %10.2f',mean(mean_pinf_bGM),mean(mean_pinf_aGM),mean(std_pinf_bGM),mean(std_pinf_aGM));
fprintf('\n');

fprintf('dy    %10.2f %10.2f %10.2f %10.2f',mean(mean_dy_bGM),mean(mean_dy_aGM),...
    mean(std_dy_bGM),mean(std_dy_aGM));
fprintf('\n');

fprintf('dc    %10.2f %10.2f %10.2f %10.2f',mean(mean_dc_bGM),mean(mean_dc_aGM),...
    mean(std_dc_bGM),mean(std_dc_aGM));
fprintf('\n');

fprintf('dinve %10.2f %10.2f %10.2f %10.2f',mean(mean_dinve_bGM),mean(mean_dinve_aGM),...
    mean(std_dinve_bGM),mean(std_dinve_aGM));
fprintf('\n');

fprintf('dw    %10.2f %10.2f %10.2f %10.2f',mean(mean_dw_bGM),mean(mean_dw_aGM),...
    mean(std_dw_bGM),mean(std_dw_aGM));
fprintf('\n');

fprintf('y                           %10.2f %10.2f',mean(std_y_bGM),mean(std_y_aGM));
fprintf('\n');

fprintf('c                           %10.2f %10.2f',mean(std_c_bGM),mean(std_c_aGM));
fprintf('\n');

fprintf('  Percentage of GM explained on average: pinf dy dc y c\n');
fprintf('%10.2f %10.2f %10.2f %10.2f %10.2f \n',GM_pinf,GM_dy,GM_dc,GM_yobs,GM_cobs);

fprintf('  Percentage of runs exceeding GM drop: pinf dy dc y c\n');
fprintf('%10.2f %10.2f %10.2f %10.2f %10.2f \n',count_GM_pinf,count_GM_dy,count_GM_dc,count_GM_yobs,count_GM_cobs);
fprintf('\n');
