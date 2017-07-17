function [N_min,N_inter,N_max]=router1_LRU_plain(produ,t_inst,ProbForSavingR1,ProbForSavingR2,N_min,N_inter,N_max)%,Sele_1,Sele_2)
global memoryR1_LRU Pop_producers Router1_hit_count FreshnessMax FreshnessMin%memoryR1_LRU memoryR1_Random;

%% Least Recently Used (LRU) Policy

% Remove data with LRU Policy. Use followings conditions for
% implementation.

% if cache empty
%     store the data
% else
%     if data exist without freshness
%         replace data with new one (CONSIDER IT AS MISS)
%     else (to store new data and remove one of the old one)
%         Choose the data which is at last location (Logic is given as below) and replace that data with the new one.
%     end
% end

% Logic for LRU removal::::::::::::::
% Sort data according to t_inst, it has been used. Remove the data which is
% having oldest t_inst for its use.

% N_min will get increamented if data is found in memoryR1 else N_inter 
% will get increamented, indicates that data is not available and need 
% retrieval from memoryR2.

% Variable discription
% memoryR1_LRU: Cache for storing data
%           column1: Producers; column2: t_stamp
% Router1_hit_count: Global variable to count Router1 hit count
% FreshnessMin/FreshnessMax: Global variable for freshness
% produ: Producer number requested.
% t_inst: time instant of request
% ProbForSavingR1: Proabbaility for saving at Router1
% ProbForSavingR2: Proabbaility for saving at Router2
% Pop_producers: Number of popular producers
% N_min,N_inter,N_max: Number of requests served by Router1, Router2 and
%                      Producers

% Remove the data which is least recently used. To implement that we had
% sorted data in decreasing order of their time_stamp of usage. Data of last
% location will get replaced by new data and again get sorted according to
% time_stamp of usage. Do not remove stale data as this is plain LRU
% implementation.

%################ TO REMOVE EXPIRED/STALE PRODUCER CONTENT ######################
% To remove expired producer content
% t_inst
% memoryR1_LRU
% indices1=find(memoryR1_LRU(:,2)>Pop_producers); % Identifying less popular users
% indices2=find(t_inst - memoryR1_LRU(indices1,3)>FreshnessMax); % Identifying stale data
% memoryR1_LRU(indices1(indices2),:)=0;
% clear indices1 indices2
% 
% indices1=find(memoryR1_LRU(:,2)<Pop_producers+1); % Identifying more popular users
% indices2=find(t_inst - memoryR1_LRU(indices1,3)>FreshnessMin); % Identifying stale data
% memoryR1_LRU(indices1(indices2),:)=0;
% clear indices1 indices2
% % memoryR1_LFU
% [~,indices1]=sort(memoryR1_LRU(:,1),'descend');
% memoryR1_LRU=memoryR1_LRU(indices1,:);
% memoryR1_LRU
%########################### END REMOVING #################################


% temp1(:,1)=memoryR1_LFU(:,1); % Considering Producers Only
index=find(memoryR1_LRU(:,2) ==produ,1,'first'); % Check for the producer
% Frequency_R1(produ,1)=Frequency_R1(produ,1)+1; % Increament the Frequency for Producer
%         index
        if ~isempty(index) % True implies producer is present
%             temp2=memoryR1_LFU(index,2);
            memoryR1_LRU(index,1)=t_inst;
%             Router1_hit_count(produ)=Router1_hit_count(produ)+1;
            if produ <= Pop_producers % check for type of producer
                if ((t_inst- memoryR1_LRU(index,3))<=FreshnessMin) % Check for freshness % && ((t_inst-memoryR1_LFU(index,3))>=0) 
%                     display('producer present At R1 with data');
                    N_min=N_min+1; % Data is present in R1
                    Router1_hit_count(produ)=Router1_hit_count(produ)+1;
                else
%                     display('producer present At R1 without data');
                    [N_inter,N_max,t_stamp]=router2_LRU_plain(produ,t_inst,ProbForSavingR2,N_inter,N_max);
                    memoryR1_LRU(index,3)=t_stamp;
                end
            else
                if ((t_inst- memoryR1_LRU(index,3))<=FreshnessMax) % Check for freshness % && ((t_inst-memoryR1_LFU(index,3))>=0) 
%                     display('producer present At R1 with data');
                    N_min=N_min+1; % Data is present in R1
                    Router1_hit_count(produ)=Router1_hit_count(produ)+1;
                else
%                     display('producer present At R1 without data');
                    [N_inter,N_max,t_stamp]=router2_LRU_plain(produ,t_inst,ProbForSavingR2,N_inter,N_max);
                    memoryR1_LRU(index,3)=t_stamp;
                end
            end
                % Sorting content in memory in decreasing order of
                % time_stamp it is being used i.e. according to column 1.
                [~,indices1]=sort(memoryR1_LRU(:,1),'descend');
                memoryR1_LRU=memoryR1_LRU(indices1,:);
%                 display('Data Present in R1')
%             end % Case where producer is present with expired data will 
                % never happen as expired data is cleared in the initial part of the code.            
        else % Case when producer is not present in CacheR1
%             display('producer not present at R1')
            [N_inter,N_max,t_stamp]=router2_LRU_plain(produ,t_inst,ProbForSavingR2,N_inter,N_max);%,Sele_2);
      % Check for empty location and index of least frequently used producer       
%             memory(:,1)=memoryR1_LFU(:,1);
%             Freshness=Frequency_R1;   
            ProbForSaving1=0;
            index2=length(memoryR1_LRU);
%             [ProbForSaving1,index2]=FindLRU(memoryR1_LRU(:,1),ProbForSavingR1);
          
    % Genrate choice variable according to probabilty ProbForSaving
            if rand()<max(ProbForSavingR1,ProbForSaving1)
                choice=1;             
            else
                choice=0;
            end
%             display('Router1 Choice');
%             choice
            if choice==1
                memoryR1_LRU(index2,:)=[t_inst,produ,t_stamp];
            end
            
            % Sorting content in memory in decreasing order of
            % time_stamp it is being used i.e. according to column 1.
            [~,indices1]=sort(memoryR1_LRU(:,1),'descend');
            memoryR1_LRU=memoryR1_LRU(indices1,:);
            
        end
%         memoryR1_LRU
        clear temp1 temp2
%         memoryR1_LRU
%         Freshness_R1
% memoryR1_LFU(size(memoryR1_LFU,1),:)=[1,produ,Frshness];
%         [~,d2]=sort(memoryR1_LFU(:,1),'descend');
%         memoryR1_LFU=memoryR1_LFU(d2,:);

end