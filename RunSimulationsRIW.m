clear all;

global n y n_ins_neg n_ins_unknown n_ins_pos n_recw_neg n_recw_unknown n_recw_pos n_rec_neg n_rec_unknown n_rec_pos c undiag_u undiag_n;

Params.TrackingCircumcision=1;
Params.TrackingWithdrawal=1;

if Params.TrackingCircumcision==1
    NumberOfParamsToOptimise=4;
elseif Params.TrackingWithdrawal==1
    NumberOfParamsToOptimise=3;
else
    NumberOfParamsToOptimise=2;
end



[prelim_data headers]=xlsread('newdataxcirc.xls');%xlsread('Data.xls');
[num_subjects, num_cols] = size(prelim_data);

index_all=0;
index_positive=0;
index_unknown_or_pos=0;
for i=1:num_subjects
    if sum(prelim_data(i,2:10))==0
    else
        index_all=index_all+1;
        data(index_all,:)=prelim_data(i,:);

        %Include those who had known HIV-positive partners
        if prelim_data(i,4)+prelim_data(i,7)+prelim_data(i,10)>0 
            index_positive=index_positive+1;
            data_discordant(index_positive,:)=prelim_data(i,:);
        end
        
        %Include those who had known HIV-positive partners or of unknown
        %status
        if sum(prelim_data(i,3:4))+sum(prelim_data(i,6:7))+sum(prelim_data(i,9:10))>0 
            index_unknown_or_pos=index_unknown_or_pos+1;
            data_discordant_unknown(index_unknown_or_pos,:)=prelim_data(i,:);
        end
    end
end
[N_all, cols]=size(data);
[N_disc, cols]=size(data_discordant);
[N_disc_unknown, cols]=size(data_discordant_unknown);


% All contacts
data_for_analysis=data;
n=N_all;

y_temp=data_for_analysis(:,11);

n_ins_neg_temp=data_for_analysis(:,4);
n_ins_unknown_temp=data_for_analysis(:,3);
n_ins_pos_temp=data_for_analysis(:,2);
n_recw_neg_temp=data_for_analysis(:,7);
n_recw_unknown_temp=data_for_analysis(:,6);
n_recw_pos_temp=data_for_analysis(:,5);
n_rec_neg_temp=data_for_analysis(:,10);
n_rec_unknown_temp=data_for_analysis(:,9);
n_rec_pos_temp=data_for_analysis(:,8);
c_temp=data_for_analysis(:,12);


NumSims=100;%10000;

uu=0;
un=0;

ResultsMatrix=zeros([NumberOfParamsToOptimise NumSims]);

tic
sim_step=1; 
TotalSims=12*NumSims;
for undiag_n=[0.005 0.01 0.015 0.02] 
un=un+1;
uu=0;
for undiag_u=[0.05 0.10 0.15]
uu=uu+1;
    

for sim_num=1:NumSims
    % Pick individuals at random for simulation (
    vec=zeros(n,1);
    for i=1:n
        vec(i)=ceil(rand()*n);
    end
    selection_index = sort(vec);
    count_freq=zeros(n,1);
    for i=1:n
        count_freq(selection_index(i))=count_freq(selection_index(i))+1;
    end
    
    row_in_data=0;
    for individual=1:n
        for freq_for_ind=1:count_freq(individual)
            row_in_data=row_in_data+1;
            y(row_in_data,1)=y_temp(individual);
            n_ins_neg(row_in_data,1)=n_ins_neg_temp(individual);
            n_ins_unknown(row_in_data,1)=n_ins_unknown_temp(individual);
            n_ins_pos(row_in_data,1)=n_ins_pos_temp(individual);
            n_recw_neg(row_in_data,1)=n_recw_neg_temp(individual);
            n_recw_unknown(row_in_data,1)=n_recw_unknown_temp(individual);
            n_recw_pos(row_in_data,1)=n_recw_pos_temp(individual);
            n_rec_neg(row_in_data,1)=n_rec_neg_temp(individual);
            n_rec_unknown(row_in_data,1)=n_rec_unknown_temp(individual);
            n_rec_pos(row_in_data,1)=n_rec_pos_temp(individual);
            c(row_in_data,1)=2-c_temp(individual);%change 2 to 0, change 1 to 1 (2 is uncircumsised, 1 is circumsised)
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %do a random walk minimisation
    max_iterations=1000000;
    
    current_position=0.1*ones(1, NumberOfParamsToOptimise);
    
    power=-2;
    p=1;
    last_move=0;
    current_height=Function_riw_beta(current_position, Params);
    
    
    while p<max_iterations && power>-6;
        
        p=p+1;
        
        %create a random amount to move by
        v=2*rand(1,NumberOfParamsToOptimise)-1;
        
        
        %change the vector size to be the right size for this step. 
        v=10^power*v;
        
        next_position=current_position+v;
        
        
        %note need to test if outside of (0, 1)
        if min(next_position)>0 && max(next_position)<1
            next_height=Function_riw_beta(next_position, Params);
            if next_height<current_height
                current_position=next_position;
                last_move=0;
                current_height=next_height;
                %dlmwrite('trace.txt',[current_position p power],'-append','newline','pc','delimiter','\t');
            else
                last_move=last_move+1;
            end
        end
        
        if last_move>100 %if its taken more than 100 tests to find a lower point, could be any number
            power=power-1; %make the vector smaller by a factor of 10
            last_move=0;%reset the last_move
        end
        
    end
    
%     hold on;
%         loglog(current_position(1),current_position(2) , '.');
%         %xlim(0, 1);
%         %ylim(0, 1);
%         drawnow;
    
    %Put the optmised result into the matrix
    ResultsMatrix(:, sim_num)=current_position;
    
    %Display where everything is up to
    if mod(sim_step,10)==0
        current_time=toc;
        fprintf('\n \n \nSimulation running: %d %d %d / %d\n', undiag_u, undiag_n, sim_step, TotalSims);
        fprintf('Time running: %f hours\n', current_time/60/60);
        time_to_finish=(TotalSims-sim_step)/(sim_step/current_time);
        fprintf('Time to finish: %f hours\n', time_to_finish/3600);
    end
    sim_step=sim_step+1;
end


ResultsStruct(uu, un).undiag_n=undiag_n;
ResultsStruct(uu, un).undiag_u=undiag_u;

ResultsStruct(uu, un).Receptive=ResultsMatrix(1, :);
ResultsStruct(uu, un).Insertive=ResultsMatrix(2, :);

i=1;
ResultsStruct(uu, un).Median(i)=median(ResultsMatrix(i, :));
ResultsStruct(uu, un).Mean(i)=mean(ResultsMatrix(i, :));
ResultsStruct(uu, un).LCI(i)=prctile(ResultsMatrix(i, :), 2.5);
ResultsStruct(uu, un).UCI(i)=prctile(ResultsMatrix(i, :), 97.5);
i=2;
ResultsStruct(uu, un).Median(i)=median(ResultsMatrix(i, :));
ResultsStruct(uu, un).Mean(i)=mean(ResultsMatrix(i, :));
ResultsStruct(uu, un).LCI(i)=prctile(ResultsMatrix(i, :), 2.5);
ResultsStruct(uu, un).UCI(i)=prctile(ResultsMatrix(i, :), 97.5);


if Params.TrackingWithdrawal==1
    ResultsStruct(uu, un).Withdrawal=ResultsMatrix(3, :);
    i=3;
    ResultsStruct(uu, un).Median(i)=median(ResultsMatrix(i, :));
    ResultsStruct(uu, un).Mean(i)=mean(ResultsMatrix(i, :));
    ResultsStruct(uu, un).LCI(i)=prctile(ResultsMatrix(i, :), 2.5);
    ResultsStruct(uu, un).UCI(i)=prctile(ResultsMatrix(i, :), 97.5);
end
if Params.TrackingCircumcision==1
    ResultsStruct(uu, un).Circumcision=ResultsMatrix(4, :);
    i=4;
    ResultsStruct(uu, un).Median(i)=median(ResultsMatrix(i, :));
    ResultsStruct(uu, un).Mean(i)=mean(ResultsMatrix(i, :));
    ResultsStruct(uu, un).LCI(i)=prctile(ResultsMatrix(i, :), 2.5);
    ResultsStruct(uu, un).UCI(i)=prctile(ResultsMatrix(i, :), 97.5);
end




end
end




