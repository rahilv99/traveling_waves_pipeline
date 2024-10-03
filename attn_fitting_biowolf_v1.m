%% 5b. Circular 2D Linear Regresssion Model Fitting for all time points and events at a 
%% specified  cluster 

function attn_fitting_biowolf_v1(ct, el, folder_name)

folder_name = fullfile(folder_name, 'Processing');

ct = str2double(ct);
el = str2double(el);
load(fullfile(folder_name, "allwave.mat"));

score = wave_info(ct).score;
relative_phase = wave_info(ct).relative_phase;
dist = wave_info(ct).dist;
dims = wave_info(ct).dims; % chan x event x time

ang_el=zeros(dims(3),dims(2));
spatial_frequency_el=zeros(dims(3),dims(2));
cl_corr_el=zeros(dims(3),dims(2));
el_dist = 25; 

for tt=1:dims(3)
    for i=1:dims(2)
        if ~all(isnan(relative_phase(i,tt,:)))

            circular = squeeze(relative_phase(i,tt,dist(el,:) < el_dist));
            linear = score(dist(el,:) < el_dist,:);

            pos_x = linear(:,1);
            pos_y = linear(:,2);

            phase = circular;

            myfun1 = @(p) -sqrt((sum(cos(phase-p(1)*pos_x-p(2)*pos_y)/length(phase)).^2 + (sum(sin(phase-(p(1)*pos_x)-p(2)*pos_y))/length(phase)).^2));

            t1=0;

            s=0;

            for theta=pi*(0:5:355)/180
                t1=t1+1;
                t2=0;
                for r=[0:.5:9 10:18 ]*pi/180
                    t2=t2+1;
                    ss(t1,t2)=myfun1([r*cos(theta),r*sin(theta)]);

                    if ss(t1,t2)<s
                        s=ss(t1,t2);

                        ThetaR=[theta,r];
                    end
                end
            end
            %     cc=s;
            %     keyboard
            sl=ThetaR(2)*[cos(ThetaR(1)) sin((ThetaR(1)))];

            ang_el(tt,i)=ThetaR(1);
            spatial_frequency_el(tt,i)=ThetaR(2);

            % calculate offset
            offs = atan2(sum(sin(phase-sl(1)*pos_x-sl(2)*pos_y)),sum(cos(phase-sl(1)*pos_x-sl(2)*pos_y))) ;

            % circular-linear correlation:
            pos_circ = mod(sl(1)*pos_x+sl(2)*pos_y+offs, 2*pi); % circular variable derived from the position
            phase_mean = mod(angle(sum(exp(1i*phase))/length(phase)),2*pi); % circular mean of the theta phase
            pos_circ_mean = mod(angle(sum(exp(1i*pos_circ))/length(phase)),2*pi); % circular mean of the circular position variable
            % calculating the correlation- this is pgd value! no adjustment for
            % number of fitted model params
            cl_corr_el(tt,i) = sum(sin(phase - phase_mean) .* sin(pos_circ - pos_circ_mean)) / sqrt( sum(sin(phase - phase_mean).^2) * sum(sin(pos_circ - pos_circ_mean).^2) );
        end
    end
end

el_data.ang_el = ang_el;
el_data.spatial_freq = spatial_frequency_el;
el_data.cl_corr = cl_corr_el;

out_dir = fullfile(folder_name, num2str(ct));

if ~exist(out_dir, 'dir') % make folder if doesn't exist already
    mkdir(out_dir);
end
save(fullfile(out_dir,num2str(el)), 'el_data', '-v7.3');

clear; close all;
pause(5);
