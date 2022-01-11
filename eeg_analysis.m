%rat eeg (burst supression ratio)
read_Intan_RHS2000_file('10032020 RAT 1_201003_123129.rhs')

data_all = ans;
data = data_all.amplifier_channels;
fs = data_all.sample_rate;
number_of_ch = data_all.num_amplifier_channels;
time = data_all.time;
const = 200*fs;

ch1 = data(1).amplifier_data;
ch2 = data(2).amplifier_data;
ch3 = data(3).amplifier_data;
ch4 = data(4).amplifier_data;



ch1_abs = abs(ch1);ch1_1 = ch1(1:fs*100*10);
ch1_2 = ch1((fs*100*10)+1:fs*100*20);
ch1_3 = ch1((fs*100*20)+1:fs*100*30);
ch1_4 = ch1((fs*100*30)+1:fs*100*40);
ch1_5 = ch1((fs*100*40)+1:fs*100*50);
ch1_6 = ch1((fs*100*50)+1:fs*100*60);

ch2_1 = ch2(1:fs*100*10);
ch2_2 = ch2((fs*100*10)+1:fs*100*20);
ch2_3 = ch2((fs*100*20)+1:fs*100*30);
ch2_4 = ch2((fs*100*30)+1:fs*100*40);
ch2_5 = ch2((fs*100*40)+1:fs*100*50);
ch2_6 = ch2((fs*100*50)+1:fs*100*60);

ch3_1 = ch3(1:fs*100*10);
ch3_2 = ch3((fs*100*10)+1:fs*100*20);
ch3_3 = ch3((fs*100*20)+1:fs*100*30);
ch3_4 = ch3((fs*100*30)+1:fs*100*40);
ch3_5 = ch3((fs*100*40)+1:fs*100*50);
ch3_6 = ch3((fs*100*50)+1:fs*100*60);

ch4_1 = ch4(1:fs*100*10);
ch4_2 = ch4((fs*100*10)+1:fs*100*20);
ch4_3 = ch4((fs*100*20)+1:fs*100*30);
ch4_4 = ch4((fs*100*30)+1:fs*100*40);
ch4_5 = ch4((fs*100*40)+1:fs*100*50);
ch4_6 = ch4((fs*100*50)+1:fs*100*60);

figure;plot(time,ch1)
subplot(4,1,1);plot(time,ch1)
subplot(4,1,2);plot(time,ch2)
subplot(4,1,3);plot(time,ch3)
subplot(4,1,4);plot(time,ch4)
%%
all_burst = zeros(4,const);
figure
for r = 1:4;
    ch1 = data(r).amplifier_data;
    subplot(4,1,r)
    [b,a] = butter(2,[3/(fs/2) 100/(fs/2)],'bandpass');
    filtered_data = filtfilt(b,a,ch1);
    % [b,a] = butter(2,[55/(fs/2) 65/(fs/2)],'stop');
    % filtered_data = filtfilt(b,a,filtered_data);
    
    trial = abs(filtered_data);
    t = time;
    i = 0;
    leng = 5000;
    while i <=const;
        
        new_ch1 = mean(trial(1+i:leng+i));
        if new_ch1>=15;
            
            trial(1+i:leng+i)=100;
            
        else
            trial(1+i:leng+i)=0;
        end
        
        if i>12001
            if new_ch1>=15 && trial(i-1)==0 && trial(i-6005)>0
                trial(1+i-leng:leng+i)=100;
            end
        end
        
        i = i+leng;
        
    end
    
    % figure;plot(time(1:fs*600),filtered_data(1:fs*600));hold on;plot(time(1:fs*600),trial(1:fs*600))
    
    trial = trial(1:const);
    %%% distance between burst
    
    start_end_points = diff(trial);
    
    [pks_start, lcs_start] = findpeaks(trial);
    
    [pks_end, lcs_end] = findpeaks(-trial);
    
    % figure;plot(trial);hold on;plot(lcs_start,pks_start,'r*');plot(lcs_end,pks_end,'k*')
    
    if lcs_start(1)>lcs_end(1);
        lcs_start = [100 lcs_start];
        pks_start = [100 pks_start];
    end
    if lcs_end(end)<lcs_start(end);
        lcs_start(end) = [];
        pks_start(end) = [];
    end
    
    silence_period = lcs_start(2:end)-lcs_end(1:end-1);
    
    burst_period = [];
    for s = 1:numel(silence_period);
        if silence_period(s)>(0.5*fs)
            burst_period = [burst_period s];
        end
    end
    
    
    lcs_start_time = lcs_start./fs;
    lcs_end_time = lcs_end./fs;
    burst_lenght = 0;
    plot(time(1:const),filtered_data(1:const));hold on;
    plot(time(1:const),trial);hold on;
    for i =1:numel(burst_period);
        all_burst(r,(lcs_start(burst_period(i)):lcs_end(burst_period(i)))) = 100;
        burst_lenght = burst_lenght+(lcs_end(burst_period(i))-lcs_start(burst_period(i)));
        plot((lcs_start(burst_period(i)):lcs_end(burst_period(i)))./fs,ones(numel(lcs_start(burst_period(i)):lcs_end(burst_period(i))),1)*100,'r','LineWidth',2);
    end
    box off;
    set(gca,'FontSize',10)
    set(gca,'FontWeight','bold')
    
%     xlim([15 30])
    
end


all_one_new = mean(all_burst);
for i =1:numel(all_one_new);
    if all_one_new(i) ~= 100;
        all_one_new(i) = 0;
    end
end

for ii=1:4;
    subplot(4,1,ii);hold on;
    area(time(1:const),all_one_new,'FaceAlpha',0.3)
end

(sum(all_one_new)/100)/numel(all_one_new)
%%
figure;plot(trial);hold on;
for i =1:numel(burst_period);
    plot((lcs_start(burst_period(i)):lcs_end(burst_period(i))),ones(numel(lcs_start(burst_period(i)):lcs_end(burst_period(i))),1),'r');
end


burst_interval = lcs_start(lcs_start(burst_period (1)):lcs_start(burst_period (1)));



