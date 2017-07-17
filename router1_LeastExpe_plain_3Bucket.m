function [N_min,N_max]=router1_LeastExpe_plain_3Bucket(produ,t_inst,ProbForSavingR1,N_min,N_max)
global memoryR1_LeastExpe Probability_producers Pop_producers Freshness_requirment Router1_hit_count count1%memoryR1_LRU memoryR1_Random;

%% Least Expected Policy

% Remove data with Least Expected Policy. Use followings conditions for
% implementation.

% if cache empty
%     store the data
% else
%     if data exist without freshness
%         replace data with new one (CONSIDER IT AS MISS)
%     else (to store new data and remove one of the old one)
%         Choose the data which is least expected in future (Logic is given as below) and replace that data with the new one.
%     end
% end

% Logic for Least expected removal::::::::::::::
% calculate E[Life-time of data]=remaining time*Probability of request in
% remaing time. Remove the data which is having least expected life-time.

% N_min will get increamented if data is found in memoryR1 else N_inter 
% will get increamented, indicates that data is not available and need 
% retrieval from memoryR2.

% Variable discription
% memoryR1_LeastExpe: Cache for storing data
%           column1: Producers; column2: t_stamp
% count1: Temprary variable for checking empty cache
%         if count1>length(cache) => Not empty
% Router1_hit_count: Global variable to count Router1 hit count
% Freshness_requirment: Global variable for freshness
% produ: Producer number requested.
% t_inst: time instant of request
% ProbForSavingR1: Proabbaility for saving at Router1
% Pop_producers: Number of type of producers
% N_min,N_max: Number of requests served by Router1 and Producers

%################ TO REMOVE EXPIRED/STALE PRODUCER CONTENT ######################
% To remove expired producer content
% t_inst
% memoryR1_LeastExpe
% indices1=find(memoryR1_LeastExpe(:,1)>Pop_producers); % Identifying less popular users
% indices2=find(t_inst - memoryR1_LeastExpe(indices1,2)>FreshnessMax); % Identifying stale data
% memoryR1_LeastExpe(indices1(indices2),:)=0;
% clear indices1 indices2
% 
% indices1=find(memoryR1_LeastExpe(:,1)<Pop_producers+1); % Identifying more popular users
% indices2=find(t_inst - memoryR1_LeastExpe(indices1,2)>FreshnessMin); % Identifying stale data
% memoryR1_LeastExpe(indices1(indices2),:)=0;
% clear indices1 indices2
% % memoryR1_LFU
% [~,indices1]=sort(memoryR1_LeastExpe(:,2),'descend');
% memoryR1_LeastExpe=memoryR1_LeastExpe(indices1,:);
% memoryR1_LeastExpe
%########################### END REMOVING #################################

% temp1(:,1)=memoryR1_LFU(:,1); % Considering Producers Only
index=find(memoryR1_LeastExpe(:,1) ==produ,1,'first'); % Check for the producer
% Frequency_R1(produ,1)=Frequency_R1(produ,1)+1; % Increament the Frequency for Producer
%         index
if ~isempty(index) % True implies producer is present
%     temp2=memoryR1_LFU(index,2);
%     memoryR1_LeastExpe(index,1)=t_inst;
%     Router1_hit_count(produ)=Router1_hit_count(produ)+1;
    temp1=sum(produ<=cumsum(Pop_producers));
    if (t_inst-memoryR1_LeastExpe(index,2))<=Freshness_requirment(temp1)
%         display('Producer present with data at R1')
        N_min=N_min+1;
        Router1_hit_count(produ)=Router1_hit_count(produ)+1;
    else % MISS HAPPENED
%         display('Producer present without data at R1')
        N_max=N_max+1;
        memoryR1_LeastExpe(index,2)=t_inst;
    end
else % Case when producer is not present in CacheR1
%     display('producer not present at R1')
    N_max=N_max+1;

    count1=count1+1;

    if count1>length(memoryR1_LeastExpe)
        [Value,index2]=FindLeastExpe_3Bucket(memoryR1_LeastExpe(:,1),memoryR1_LeastExpe(:,2),t_inst);
        if (Probability_producers(produ)*Freshness_requirment(sum(produ<=cumsum(Pop_producers))))<Value
            index2=0;
        end
    else
        index2=count1;
    end


    ProbForSaving1=0;
% Genrate choice variable according to probabilty ProbForSaving
    if rand()<max(ProbForSavingR1,ProbForSaving1)
        choice=1;             
    else
        choice=0;
    end
%             index2
%             display('Router1 Choice');
%             choice
    if choice==1 && index2~=0
        memoryR1_LeastExpe(index2,:)=[produ,t_inst];
    end

end
% [~,indices1]=sort(memoryR1_LeastExpe(:,2),'descend');
% memoryR1_LeastExpe=memoryR1_LeastExpe(indices1,:);

% memoryR1_LeastExpe
clear temp1 temp2

end