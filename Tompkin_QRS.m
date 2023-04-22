close all;
clear all;
clc;

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end
SAVE = 0;
%% 

fs = 200;
file = '7.txt';
dat = load(file);
data = dat(1:2000); % work with 2000 samples
t = (1 : 2000)/fs;
figure(10);
    plot(t, data - mean(data), 'k');
        title([file ' ECG data']);
        grid on;
        ylabel('$y(t)$');
        xlabel('$t$ [s]');
QRS_width = (2.25 - 2.18)*1.05; %[s] From visual inspection
N = ceil(QRS_width*fs);
if(SAVE)
    f = figure(10);
        f.Name = [file ' ECG data'];
        saveas(f,['.\izvestaj\slike\' f.Name '.eps'],'epsc');

end

%% 
alpha = 0.99;
fraction = 0.13;
beat = my_Pan_Tompkin_QRS(file, alpha, N, fraction);

if(SAVE)
    f = figure(1);
        f.Name = ['res_filt_' file(1)];
        saveas(f,['.\izvestaj\slike\' f.Name '.eps'],'epsc');
    f = figure(2);
        f.Name = ['peak_loc_' file(1)];
        saveas(f,['.\izvestaj\slike\' f.Name '.eps'],'epsc');
end
%%
%%
N = 50;

alpha = 0.99;
fraction = 0.12;
beat = my_Pan_Tompkin_QRS(file, alpha, N, fraction);

if(SAVE)
    f = figure(1);
        f.Name = ['res_filt_' file(1)];
        saveas(f,['.\izvestaj\slike\' [f.Name num2str(N)] '.eps'],'epsc');
    f = figure(2);
        f.Name = ['peak_loc_' file(1)];
        saveas(f,['.\izvestaj\slike\' [f.Name num2str(N)] '.eps'],'epsc');
end
%%
fs = 200;
file = '8.txt';
dat = load(file);
data = dat(1:2000); % work with 2000 samples
t = (1 : 2000)/fs;
figure(10);
    plot(t, data - mean(data), 'k');
        title([file ' ECG data']);
        grid on;
        ylabel('$y(t)$');
        xlabel('$t$ [s]');
QRS_width = (2.695 - 2.56)*1.05; %[s] From visual inspection
N = ceil(QRS_width*fs);
if(SAVE)
    f = figure(10);
        f.Name = [file ' ECG data'];
        saveas(f,['.\izvestaj\slike\' f.Name '.eps'],'epsc');

end
%%
alpha = 0.99;
fraction = 0.13;
beat = my_Pan_Tompkin_QRS(file, alpha, N, fraction);

if(SAVE)
    f = figure(1);
        f.Name = ['res_filt_' file(1)];
        saveas(f,['.\izvestaj\slike\' f.Name '.eps'],'epsc');
    f = figure(2);
        f.Name = ['peak_loc_' file(1)];
        saveas(f,['.\izvestaj\slike\' f.Name '.eps'],'epsc');
end
%%
N = 50;

alpha = 0.99;
fraction = 0.13;
beat = my_Pan_Tompkin_QRS(file, alpha, N, fraction);

if(SAVE)
    f = figure(1);
        f.Name = ['res_filt_' file(1)];
        saveas(f,['.\izvestaj\slike\' [f.Name num2str(N)] '.eps'],'epsc');
    f = figure(2);
        f.Name = ['peak_loc_' file(1)];
        saveas(f,['.\izvestaj\slike\' [f.Name num2str(N)] '.eps'],'epsc');
end
%%

%Low-Pass
B1 = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
A1 = [1 -2 1];
[Hlp, w] = freqz(B1, A1);

% HP filter (check Equation 2.3 in the Errata)
%
    B2 = zeros(1, 32);
    B2(1) = -1/32; B2(16) = 1; B2(17) = -1; B2(32) = 1/32;
    A2 = [1 -1];
[Hhp, w] = freqz(B2, A2, w);
w = w/pi*0.5*fs;


figure(100);

    subplot(3, 1, 1);
    semilogx(w, 20*log(Hlp), 'k');
    hold on;
        title('Low-pass filter, frequency response');
        ylabel('$\|H(f)\|$');
        grid on;
        xlim([w(1) w(end)]);
        
    subplot(3, 1, 2);
    semilogx(w, 20*log10(Hhp), 'k');
    hold on;
        title('High-pass filter, frequency response');
        ylabel('$\|H(f)\|$');
        grid on;
        xlim([w(1) w(end)]);
        
    subplot(3, 1, 3);
    semilogx(w, 30*log((Hhp.*Hlp)), 'k');
    hold on;
        title('Band-pass filter, frequency response');
        ylabel('$\|H(f)\|$');
        xlabel('$f$ [Hz]');
        grid on;
        xlim([w(1) w(end)]);  
        
 if(SAVE)
    f = figure(100);
        f.Name = 'filter';
        saveas(f,['.\izvestaj\slike\' f.Name '.eps'],'epsc');
end