function [N_inter,N_max,t_stamp]=router2_LRU_plain(produ,t_inst,ProbForSavingR2,N_inter,N_max)%,Sele_2)
global memoryR2_LRU Pop_producers Router2_hit_count FreshnessMin FreshnessMax%memoryR2_Random;
% global LocInFile;

% CacheSize=length(memoryR2_LFU);
% filename=sprintf('TestdataCacheSize%d.xlsx',CacheSize);

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

% N_inter will get increamented if data is found in memoryR2 else N_max 
% will get increamented, indicates that data is not available and need 
% retrieval from producer.

% Variable discription
% memoryR2_LRU: Cache for storing data
%           column1: Producers; column2: t_stamp
% Router2_hit_count: Global variable to count Router1 hit count
% FreshnessMin/FreshnessMax: Global variable for freshness
% produ: Producer number requested.
% t_inst: time instant of request
% ProbForSavingR2: Proabbaility for saving at Router2
% Pop_producers: Number of popular producers
% N_inter,N_max: Number of requests served by Router2 and
%                      Producers

% Remove the data which is least recently used. To implement that we had
% sorted data in decreasing order of their time_stamp of usage. Data of last
% location will get replaced by new data and again get sorted according to
% time_stamp of usage. Do not remove stale data as this is plain LRU
% implementation.

%################ TO REMOVE EXPIRED PRODUCER CONTENT ######################
    % To remove expired producer content
% t_inst   
% memoryR2_LRU
% indices1=find(memoryR2_LRU(:,2)>Pop_producers); % Identifying less popular users
% indices2=find(t_inst - memoryR2_LRU(indices1,3)>FreshnessMax); % Identifying stale data
% memoryR2_LRU(indices1(indices2),:)=0;
% clear indices1 indices2
% 
% indices1=find(memoryR2_LRU(:,2)<Pop_producers+1); % Identifying more popular users
% indices2=find(t_inst - memoryR2_LRU(indices1,3)>FreshnessMin); % Identifying stale data
% memoryR2_LRU(indices1(indices2),:)=0;
% clear indices1 indices2
% % memoryR2_LFU
% [~,indices1]=sort(memoryR2_LRU(:,1),'descend');
% memoryR2_LRU=memoryR2_LRU(indices1,:);
% memoryR2_LRU
%########################### END REMOVING #################################

% temp1(:,1)=memoryR2_LFU(:,1);
index=find(memoryR2_LRU(:,2)==produ,1,'first');% Check for producer data
% Frequency_R2(produ,1)=Frequency_R2(produ,1)+1;
%         index
if ~isempty(index)
%     temp2=memoryR2_LRU(index,2);
    N_inter=N_inter+1;
    memoryR2_LRU(index,1)=t_inst;
%     Router2_hit_count(produ)=Router2_hit_count(produ)+1;
    if produ<=Pop_producers % check for type of producer
        if (t_inst-memoryR2_LRU(index,3))<=FreshnessMin
%             display('Producer present at R2 with data');
            t_stamp=memoryR2_LRU(index,3);
            Router2_hit_count(produ)=Router2_hit_count(produ)+1;
        else
%             display('Producer present at R2 without data');
%             memoryR2_LRU(index,1)=t_inst;
            memoryR2_LRU(index,3)=t_inst;
            t_stamp=t_inst;
        end
    else
        if (t_inst-memoryR2_LRU(index,3))<=FreshnessMax
%             display('Producer present at R2 with data');
            t_stamp=memoryR2_LRU(index,3);
            Router2_hit_count(produ)=Router2_hit_count(produ)+1;
        else
%             display('Producer present at R2 without data');
            memoryR2_LRU(index,1)=t_inst;
            memoryR2_LRU(index,3)=t_inst;
            t_stamp=t_inst;
        end
    end
        % Sorting content in memory in decreasing order of
        % time_stamp it is being used i.e. according to column 1.
        [~,indices1]=sort(memoryR2_LRU(:,1),'descend');
        memoryR2_LRU=memoryR2_LRU(indices1,:);
        
%                 display('producer present At R2 with data')
        clear temp2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THIS PART OF THE CODE BULID FOR CHECK%%%%%%%%%%%               
%                 temp1={'t_inst', t_inst; 'Producer',produ;'Freshness',Frshness};
%                 temp2=sprintf('A%d',LocInFile);
%                 xlswrite(filename,temp1,1,temp2);
%                 LocInFile=LocInFile+4;
%                 clear temp1 temp2
%                 
%                 temp1={'memoryR1_LFU'};
%                 temp2=sprintf('A%d',LocInFile);
%                 xlswrite(filename,temp1,1,temp2);
%                 temp1={'memoryR2_LFU'};
%                 temp2=sprintf('D%d',LocInFile);
%                 xlswrite(filename,temp1,1,temp2);
%                 LocInFile=LocInFile+2;
%                 clear temp1 temp2
%                 
%                 temp1={'Producers','T_Stamp'};
%                 temp2=sprintf('A%d',LocInFile);
%                 xlswrite(filename,temp1,1,temp2);
%                 temp2=sprintf('D%d',LocInFile);
%                 xlswrite(filename,temp1,1,temp2);
%                 LocInFile=LocInFile+1;
%                 clear temp1 temp2
%                                 
%                 temp2=sprintf('A%d',LocInFile);
%                 xlswrite(filename,memoryR1_LFU,1,temp2);
%                 temp2=sprintf('D%d',LocInFile);
%                 xlswrite(filename,memoryR2_LFU,1,temp2);
%                 
%                 LocInFile=LocInFile+CacheSize+2;
% %                 
%                 f = warndlg('This is a warning.', 'A Warning Dialog');
%                 disp('This prints immediately');
%                 drawnow     % Necessary to print the message
%                 waitfor(f);
%                 disp('This prints after you close the warning dialog');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%S
%%% THIS WILL NEVER HAPPEN AS WE ARE CLEARING DATA OF PRODUCERS NOT
%%% SATISFYING FRESHNESS CONDITIONS. AND THIS ELSE PART CORRESPONDS TO
%%% PRODUCERS WHO ARE IN CACHE BUT HAVING EXPIRED DATA.
%             else
%                 display('producer present at R2 without data')
%                 N_max=N_max+1;
%                 memoryR2_LFU(index,2)=t_inst;
%                 t_stamp=t_inst;
%     end
else
%     display('producer not present at R2')
    N_max=N_max+1;
    t_stamp=t_inst;
    % Check for empty location and index of least frequently used producer 
%             memory(:,1)=memoryR2_LFU(:,1);
%             Freshness=Frequency_R2; 
    ProbForSaving1=0;
    index2=length(memoryR2_LRU);
%             [ProbForSaving1,index2]=FindLRU(memoryR2_LRU(:,1),ProbForSavingR2);
%             index2
    if rand()<max(ProbForSavingR2,ProbForSaving1)
        choice=1;
    else
        choice=0;
    end
%             display('Router2 Choice');
%             choice
    if choice==1                
        memoryR2_LRU(index2,:)=[t_inst,produ,t_stamp];
    end

    % Sorting content in memory in decreasing order of
    % time_stamp it is being used i.e. according to column 1.
    [~,indices1]=sort(memoryR2_LRU(:,1),'descend');
    memoryR2_LRU=memoryR2_LRU(indices1,:);
end
%     memoryR2_LRU    
    clear temp1 temp2
        
%         memoryR2_LRU
%         Freshness_R2
%         [~,d2]=sort(memoryR2_LFU(:,1),'descend');
%         memoryR2_LFU=memoryR2_LFU(d2,:);

end