%Generating Toy Data
roi_data = rand(20,5,50,2);
behavior = randi([1 10],20);

%-----------Temporal ISC-----------%

%Step one: generating mean TS for each subject & ROI
mean_cell = mean_sub_roi(roi_data);


%----Pairwise Temporal ISC----%
pairwise_temporal_ISC = zeros(20,20,5);
%mean_cell = mean_cell';
for i = 1:5
    for j = 1:20 
        cell_roi = mean_cell(:,i);
        A = [cell_roi{:}]; % for each ROI, the 50 time points of 20 subs
        C = corr(A);
        pairwise_temporal_ISC(:,:,i)= C;
    end
end


%----Leave-one-out Temporal ISC----%
loo_temporal_ISC = zeros(20,5);
index_sub = (1:20);
mat_rest_sub = zeros(50,19,5);
loo_corr = zeros(1,5);
for i = 1:20
    
    index = setdiff(index_sub, i);
    cell_current_sub = mean_cell(i,:); %current subject, 5 ROIs, 50 time points
    mat_current_sub = [cell_sub{:}]; %50 by 5
    
    cell_rest_sub = mean_cell(index,:); %remaining subjects
    
    for j = 1:5
        each_roi = cell_rest_sub(:,j);
        mat_each_roi = [each_roi{:}];
        
        mat_rest_sub(:,:,j) = mat_each_roi; % 50 by 19 by 5; 5 ROIs combined
        rest_sub_mean = squeeze(mean(mat_rest_sub,2); %taking the avg of 19 remaining subs, 50 by 5
        
        
        
        loo_corr(1,j)= corr(mat_current_sub(:,j),rest_sub_mean(:,j));
        loo_temporal_ISC(i,:) = loo_corr;

    end
end


%----Dynamic ISC----%
loo_dynamic_ISC = zeros(20,5,5);
window = 10;
step = 10;
index_sub = (1:20);
current_win_holder = zeros(10,5);
rest_win_holder = zeros(10,5,19);
for i = 1:5
    for j = 1:20 
        %cell_roi = mean_cell(:,i);
        %A = [cell_roi{:}];
        index = setdiff(index_sub, j);
        current_sub = A(:,i);
        rest_sub = A (:,index);
        %Apply sliding window
        for k = 1:5
            current_win = current_sub((10*k-9):10*k,1);
            current_win_holder(:,k) =current_win; %10 by 5
            %current_win_holder = current_win_holder'; %5 by 10

            for h = 1:19
                rest_c = rest_sub(:,h);
                rest_win = rest_c((10*k-9):10*k,1);
                rest_win_holder(:,k,h)= rest_win; % 10 by 5 by 19
                %rest_win_holder = reshape(rest_win_holder,[5,10,19]);
                rest_win_mean = mean(rest_win_holder,3);
            end
         test = corr(current_win_holder,rest_win_mean);
         loo_dynamic_ISC(j,:,:) = test;
        end
    end
end
%Plot
s = size(loo_dynamic_ISC);
X_max = s(2);
X = (1:X_max);
Y = squeeze(loo_dynamic_ISC(1,1,:));
fig = plot(X,Y);
saveas(fig,'loo_dynamic_ISC_11.png');


%-----------Spatial ISC-----------%
loo_spatial_ISC = zeros(20,1);
index_sub = (1:20);
mean_mat = zeros(20,5);

for i = 1:20 
    index = setdiff(index_sub, i);
    
    for j = 1:5
        cell = mean_cell(i,j);
        mean_mat(i,j) = mean(cell{1,1});
    end
end


for i = 1:20 
    index = setdiff(index_sub, i);
    
    for j = 1:5
        current_S = mean_mat(i,:);
        current_S = current_S';
        rest_S = mean_mat(index,:);
        rest_S_mean = mean(rest_S,1);
        rest_S_mean =rest_S_mean';
        loo_spatial = corr(current_S,rest_S_mean);
        loo_spatial_ISC(i,:) = loo_spatial;
    end
end


%-----------intra-subject correlation-----------%
intrasubject_temporal_ISC = zeros(20,5);
for i = 1:20
    for j = 1:5
        sub_roi_S1 = squeeze(roi_data(i,j,:,1));
        sub_roi_S2 = squeeze(roi_data(i,j,:,2));
        corr12 = corr(sub_roi_S1,sub_roi_S2);
        intrasubject_temporal_ISC(i,j) = corr12;
    end
end


%-----------ISFC-----------%
loo_ISFC = zeros(20,5,5);
index_sub = (1:20);
mat_rest_sub = zeros(50,19,5);
for i = 1:20
    
    index = setdiff(index_sub, i);
    cell_current_sub = mean_cell(i,:); %current subject, 5 ROIs, 50 time points
    mat_current_sub = [cell_sub{:}]; %50 by 5
    current_FC = corrcoef(mat_current_sub);
    
    cell_rest_sub = mean_cell(index,:); %remaining subjects
    
    for j = 1:5
        each_roi = cell_rest_sub(:,j);
        mat_each_roi = [each_roi{:}];
        
        mat_rest_sub(:,:,j) = mat_each_roi; % 50 by 19 by 5; 5 ROIs combined
        rest_sub_mean = squeeze(mean(mat_rest_sub,2)); %taking the avg of 19 remaining subs, 50 by 5
        rest_FC = corrcoef(rest_sub_mean);
        
        ISFC = corr(current_FC,rest_FC);
        loo_ISFC(i,:,:) = ISFC;
    end
end


%-----------Pairwise with Behavior-----------%
Behavior_mat = corr(behavior);
BB_corr = zeros(20,20,5);
for i = 1:5
    roi = pairwise_temporal_ISC(:,:,i);
    BB_corr(:,:,i) = corr(roi,Behavior_mat);
end
%roi_5 = BB_corr(:,:,5)
%imshow(roi_5)
%The behavioral matrix (behavior) consists of an integer from 1 to 10 for each
%subject. This integer could represent the working memoery capacity K for each
%subject (the maximum number of items one can remember). Among all paris of subjects, 
%the higher the correlation in the behavioral correlation matrix (Behavior_mat), 
%the more similar the number of items the subjects can remember. Among all 5 ROIs, 
%the higher the overall brain-behavior correaltion is, the more this ROI represents K similarity. 
%ROI_5 shows the highest overall coorelation value.


%----Helpful Function----%
function mean_cell = mean_sub_roi(roi_data) 
          
         mean_cell = cell(20,5);
         sz = size(roi_data);
         sub = sz(1);
         r = sz(2);
         t = sz(3);
         ses = sz(4);
         session_all = zeros(t,1,ses);
         for i = 1:sub
             for j = 1:r
                 for k = 1:ses
                     corr_k = squeeze(roi_data(i,j,:,k));
                     session_all(:,:,k) = corr_k;
                     session_mean = mean(session_all, 3);
                     mean_cell(i,j) = {session_mean};
                 end
             end    
         end    
end

