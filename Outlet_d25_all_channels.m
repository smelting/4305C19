%% Receiver Loop For Channel 1 "ON" and "OFF"
clear

press_only = 0;

ds = 25;

on_ref1 = load('on_var1.mat');
on_ref1 = on_ref1.received_dup;
on_ref1 = [on_ref1' 0 0 0];
on_ref1 = downsample(on_ref1', ds);

off_ref1 = load('off_var1.mat');
off_ref1 = off_ref1.received_dup;
off_ref1 = [off_ref1' 0 0 0];
off_ref1 = downsample(off_ref1', ds);

on_ref2 = load('on_var2.mat');
on_ref2 = on_ref2.received_dup;
on_ref2 = [on_ref2' 0 0 0];
on_ref2 = downsample(on_ref2', ds);

off_ref2 = load('off_var2.mat');
off_ref2 = off_ref2.received_dup;
off_ref2 = [off_ref2' 0 0 0];
off_ref2 = downsample(off_ref2', ds);

on_ref3 = load('on_var3.mat');
on_ref3 = on_ref3.received_dup;
on_ref3 = [on_ref3' 0 0 0];
on_ref3 = downsample(on_ref3', ds);

off_ref3 = load('off_var3.mat');
off_ref3 = off_ref3.received_dup;
off_ref3 = [off_ref3' 0 0 0];
off_ref3 = downsample(off_ref3', ds);

x=1; %counter used for averaging

%rx = sdrrx('Pluto','CenterFrequency', 433.975e6, 'SamplesPerFrame', 56000, 'BasebandSampleRate', 1e6);
rx = sdrrx('Pluto','CenterFrequency', 433.92e6, 'SamplesPerFrame', 56000, 'BasebandSampleRate', 1e6);

%% Loop to Find Current Transmitted Signal
while 1
%rx = sdrrx('Pluto','CenterFrequency', 433.975e6, 'SamplesPerFrame', 56000, 'BasebandSampleRate', 1e6);
received_dup = rx();
received_dup1 = received_dup;

%Threshold value to differentiate between a pseudo "1" and "0"
threshold = 15;
received_dup(abs(real(received_dup))<threshold)=0;
received_dup(abs(real(received_dup))>=threshold)=1;
received_dup_before_avg = received_dup;
%% test
trigger1 = 0;
for i=51:length(received_dup)-50
    
if and(received_dup(i)==1,trigger1 < 1)
    trigger1 = i;
    i = i+50;
end
    
if and(sum(received_dup(i-50:i))==0,trigger1 > 1)
    received_dup(trigger1:i-50)=1;
    trigger1 = 0;   
end
end
dup_Shift = received_dup;
received_dup = downsample(received_dup, ds);

max_sum_on1=0;
max_sum_off1=0;
max_sum_on2=0;
max_sum_off2=0;
max_sum_on3=0;
max_sum_off3=0;

bound = 18650/ds-1;

for i=1:length(received_dup)-bound

corr_on1 = sum(received_dup(i:i+bound)==on_ref1);
corr_off1 = sum(received_dup(i:i+bound)==off_ref1);
corr_on2 = sum(received_dup(i:i+bound)==on_ref2);
corr_off2 = sum(received_dup(i:i+bound)==off_ref2);
corr_on3 = sum(received_dup(i:i+bound)==on_ref3);
corr_off3 = sum(received_dup(i:i+bound)==off_ref3);

if corr_on1>max_sum_on1
max_sum_on1 = corr_on1;
end

if corr_off1>max_sum_off1
max_sum_off1 = corr_off1;
end

if corr_on2>max_sum_on2
max_sum_on2 = corr_on2;
end

if corr_off2>max_sum_off2
max_sum_off2 = corr_off2;
end

if corr_on3>max_sum_on3
max_sum_on3 = corr_on3;
end

if corr_off3>max_sum_off3
max_sum_off3 = corr_off3;
end

end

final_on_corr1 = max_sum_on1;
final_off_corr1 = max_sum_off1;
final_on_corr2 = max_sum_on2;
final_off_corr2 = max_sum_off2;
final_on_corr3 = max_sum_on3;
final_off_corr3 = max_sum_off3;

max_val_arr=[];

%corr_bound = round(0.724 * (bound+1));
corr_bound = round(0.67 * (bound+1));
final_corrs = [final_on_corr1 final_off_corr1 final_on_corr2 final_off_corr2 final_on_corr3 final_off_corr3];
if sum(final_corrs>corr_bound)==6
    max_val = find(final_corrs == max(final_corrs));
    if press_only == 1
        max_val(1)=7;
    end
    max_val_arr = [max_val_arr max_val(1)];
    x = x+1;
    if x==9
        max_val_in = mode(max_val_arr);
        max_val_arr=[];
        x=1;    
    switch max_val_in
        case 1
            result = 'C1 ON'
        case 2
            result = 'C1 OFF'
        case 3
            result = 'C2 ON'
        case 4
            result = 'C2 OFF'
        case 5
            result = 'C3 ON'
        case 6
            result = 'C3 OFF'
        otherwise 
            result = 'Pressed'
            
    end
    end
else
    result = 'No Signal'
end

end
%% plot original
plt_int = 1:length(real(received_dup1))-1;
plot(plt_int, real(received_dup1(plt_int)));
ylim([-500 500]);

%% plot non averaged binary conversion (demodulation)
plt_int = 1:length(real(received_dup_before_avg))-1;
plot(plt_int, real(received_dup_before_avg(plt_int)));
ylim([-2 2]);

%% plot averaged upsampled (decode)
plt_int = 1:length(real(dup_Shift))-1;
plot(plt_int, real(dup_Shift(plt_int)));
ylim([-2 2]);
clear plot
%% plot averaged downsampled (downsampled decoded)
plt_int = 1:length(received_dup)-1;
plot(plt_int, received_dup(plt_int));
ylim([-1 2]);

%% plot ref
refs = [on_ref1 off_ref1 on_ref2 off_ref2 on_ref3 off_ref3];
ref = refs(:,max_val_in);
plt_int = 1:length(ref);
plot(plt_int, ref);
ylim([-0.5 1.5]);