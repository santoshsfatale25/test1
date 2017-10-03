%% 
% 3 Bucket model and Zipf distribution with three types of users and
% Freshness requirement as follows:
% Freshness=[F_a F_b F_c] where F_b=F_a*10^2 or some higher value.
% Setting is as follows:
% 
% 3 Bucket Uniform Distribution: 
% Prob_a, Prob_b=0.99*Prob_a, Prob_c=1-Prob_a-Prob_b
% Freshness requirement as above
% Number of Producers in each bucket are N_a=C,N_b=C,N_c=N-N_a-N_b where C
% is cache size and N is total number of producers.
% 
% Zipf Distribution:
% Use moderate Zipf parameter beta for probability distribution and use
% remaining setting as 3 Bucket Uniform distribution.
%% Probabilistic Save Implementation
clear all;
close all;
% clc;
%%
NumberOfRequests=10^5;
NumberOfIterations=10^2;
Producers=4*10^2; % Number of Producers
global Pop_producers

global Freshness_requirment
const=20;
% F_a=5;
% F_b=const*F_a;
% F_c=F_a;
% Freshness_requirment=[F_a F_b F_c];

global Router1_hit_count

ProbForSavingVectorR1=1;%0.2:0.2:1.0;%1.0;
CacheSize=10:5:40;

Prob_a=0.4;%0.25:0.05:0.45;
beta=0.8;%0.5:0.3:1.7;
hit_rate_total_Sim_Uni_LeastExpe=zeros(NumberOfIterations,length(CacheSize));
hit_rate_total_Sim_Uni_SMP=zeros(NumberOfIterations,length(CacheSize));
hit_rate_total_Sim_Uni_RAND=zeros(NumberOfIterations,length(CacheSize));
hit_rate_total_Sim_Uni_LRU=zeros(NumberOfIterations,length(CacheSize));
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

% Characteristic_time_Uni=zeros(Producers,length(CacheSize));
% Characteristic_time_Zipf=zeros(Producers,length(CacheSize));

global memoryR1_LeastExpe memoryR1_LRU memoryR1_SMP memoryR1_RAND Probability_producers

global count1 % Checks cache is empty or not.


% Exponential inter-arrival time
time=cumsum(exprnd(1,NumberOfRequests,1));
%% ######################################### 3 class Uniform Distribution ################################################

tic;
Prob_b=0.99*Prob_a;
Prob_c=ones(1,length(Prob_a))-Prob_a-Prob_b;
ProducersProbability_Uni=zeros(Producers,length(CacheSize));
% Freshness_Uni=zeros(Producers,length(CacheSize));
for ii=1:length(CacheSize)
    N_a=CacheSize(ii);
    N_b=CacheSize(ii);
    N_c=Producers-N_a-N_b;
    ProducersProbability_Uni(:,ii)=[repmat(Prob_a./N_a,N_a,1);repmat(Prob_b./N_b,N_b,1);repmat(Prob_c./N_c,N_c,1)];
%     Freshness_Uni(:,cache)=[repmat(F_a,N_a,1);repmat(F_b,N_b,1);repmat(F_c,N_c,1)];
end

Freshness_Uni=zeros(Producers,length(CacheSize));
for ii=1:length(CacheSize)
    Freshness_Uni(:,ii)=const*ProducersProbability_Uni(1,ii)./ProducersProbability_Uni(:,ii);
end

for kk=1:NumberOfIterations
    producersRequest_Uni=zeros(NumberOfRequests,length(CacheSize));
    for ii=1:length(CacheSize)
        producersRequest_Uni(:,ii)=datasample(1:Producers,NumberOfRequests,'Weights',ProducersProbability_Uni(:,ii));
    end

    Freshness_const=Freshness_Uni.*ProducersProbability_Uni;

    %% Least Expected Variables :::::::::::::::::::::::::::::::::::::::::::::::
    R1_hit_count_Uni_LeastExpe=zeros(Producers,length(CacheSize));

    R1_hit_count_Zipf_LeastExpe=zeros(Producers,length(CacheSize));

    N_min_3class_LeastExpe=zeros(length(CacheSize),1);
    N_max_3class_LeastExpe=zeros(length(CacheSize),1);

    N_min_Zipf_LeastExpe=zeros(length(CacheSize),1);
    N_max_Zipf_LeastExpe=zeros(length(CacheSize),1);


    %% 3 class Uniform Distribution Least Expected

    for nn=1:length(Prob_a)
        N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
        N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
        for cache=1:length(CacheSize)
            Probability_producers(1,:)=ProducersProbability_Uni(:,cache);
            N_a=CacheSize(cache);
            N_b=CacheSize(cache);
            N_c=Producers-N_a-N_b;
            Pop_producers=[N_a N_b N_c];
            Freshness_requirment=Freshness_Uni(:,cache);

            for jj=1:length(ProbForSavingVectorR1)
                ProbForSavingR1=ProbForSavingVectorR1(jj);
    % memoryR1_LeastExpe and memoryR2_LeastExpe have following structure.
    % First Column: Producer number
    % Second Column: t_inst at which it is being fetched from producer

                memoryR1_LeastExpe=zeros(CacheSize(cache),2);

                Router1_hit_count=zeros(Producers,1);

                count1=0;

                message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f Probability_a=%f'...
                                ,CacheSize(cache),ProbForSavingR1,Prob_a(nn));
                h=msgbox(message);
                clear message

                N_min_temp=0;
                N_max_temp=0;
                display('3 class Uniform Distribution');
                for ii=1:length(time)
    %                 display(t_inst);
                    produ=producersRequest_Uni(ii,cache);
                    t_inst=time(ii);
    %                 Frshness=Freshness(t_inst,1);
                    [N_min_temp,N_max_temp]=router1_LeastExpe_plain_3class(produ,t_inst,...
                                                                ProbForSavingR1,N_min_temp,N_max_temp);

                end
                N_min(jj,cache)=N_min_temp;
                N_max(jj,cache)=N_max_temp;

                delete(h);
                clear('h');
            end        

            R1_hit_count_Uni_LeastExpe(:,cache)=Router1_hit_count;
        end
        N_min_3class_LeastExpe(:,nn)=N_min;
        N_max_3class_LeastExpe(:,nn)=N_max;
    end
    toc

    % delete(h1);
    % clear('h1');
    display('Done Uniform');


    requests_Uni=zeros(Producers,length(CacheSize));
    for cache=1:length(CacheSize)
        for ii=1:NumberOfRequests
            requests_Uni(producersRequest_Uni(ii,cache),cache)=requests_Uni(producersRequest_Uni(ii,cache),cache)+1;
        end
    end


    % requests_Uni(1:3,1);
    % R1_hit_count_Uni(1:3,1)
    % temp1(:,1)=R1_hit_count(1:3,1);
    % temp2(:,1)=requests(1:3,1);
    % display('Hit rate Simulation for  2 class Uniform Distribution')
    % temp1./temp2
    % temp2=repmat(requests_Uni,1,length(CacheSize));
    hit_rate_Simul_Uni_LeastExpe=R1_hit_count_Uni_LeastExpe./requests_Uni;
    clear temp1 temp2
    hit_rate_total_Sim_Uni_LeastExpe(kk,:)=sum(R1_hit_count_Uni_LeastExpe)/NumberOfRequests;
end

%% SMP Variables ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
R1_hit_count_Uni_SMP=zeros(Producers,length(CacheSize));

R1_hit_count_Zipf_SMP=zeros(Producers,length(CacheSize));

N_min_3class_SMP=zeros(length(CacheSize),1);
N_max_3class_SMP=zeros(length(CacheSize),1);

N_min_Zipf_SMP=zeros(length(CacheSize),1);
N_max_Zipf_SMP=zeros(length(CacheSize),1);

%% 3 class Uniform Distribution SMP (Least Frequently Used)

for nn=1:length(Prob_a)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    for cache=1:length(CacheSize)
        Probability_producers(1,:)=ProducersProbability_Uni(:,cache);
        N_a=CacheSize(cache);
        N_b=CacheSize(cache);
        N_c=Producers-N_a-N_b;
        Pop_producers=[N_a N_b N_c];
        Freshness_requirment=Freshness_Uni(:,cache);
        
        for jj=1:length(ProbForSavingVectorR1)
            ProbForSavingR1=ProbForSavingVectorR1(jj);
% memoryR1_LeastExpe and memoryR2_LeastExpe have following structure.
% First Column: Producer number
% Second Column: t_inst at which it is being fetched from producer

            memoryR1_SMP=zeros(CacheSize(cache),2);

            Router1_hit_count=zeros(Producers,1);

            count1=0;
            count2=0;

            message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f Probability_a=%f'...
                            ,CacheSize(cache),ProbForSavingR1,Prob_a(nn));
            h=msgbox(message);
            clear message

            N_min_temp=0;
            N_max_temp=0;
            display('3 class Uniform Distribution');
            for ii=1:length(time)
%                 display(t_inst);
                produ=producersRequest_Uni(ii,cache);
                t_inst=time(ii);
%                 Frshness=Freshness(t_inst,1);
                [N_min_temp,N_max_temp]=router1_SMP_3class(produ,t_inst,...
                                                            ProbForSavingR1,N_min_temp,N_max_temp);

            end
            N_min(jj,cache)=N_min_temp;
            N_max(jj,cache)=N_max_temp;

            delete(h);
            clear('h');
        end        

        R1_hit_count_Uni_SMP(:,cache)=Router1_hit_count;
    end
    N_min_3class_SMP(:,nn)=N_min;
    N_max_3class_SMP(:,nn)=N_max;
end
toc

% delete(h1);
% clear('h1');
display('Done Uniform');

% requests_Uni(1:3,1);
% R1_hit_count_Uni(1:3,1)
% temp1(:,1)=R1_hit_count(1:3,1);
% temp2(:,1)=requests(1:3,1);
% display('Hit rate Simulation for  2 class Uniform Distribution')
% temp1./temp2
% temp2=repmat(requests_Uni,1,length(CacheSize));
hit_rate_Simul_Uni_SMP=R1_hit_count_Uni_SMP./requests_Uni;
clear temp1 temp2
hit_rate_total_Sim_Uni_SMP=sum(R1_hit_count_Uni_SMP)/NumberOfRequests;


%% RAND Variables ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
R1_hit_count_Uni_RAND=zeros(Producers,length(CacheSize));

R1_hit_count_Zipf_RAND=zeros(Producers,length(CacheSize));

N_min_3class_RAND=zeros(length(CacheSize),1);
N_max_3class_RAND=zeros(length(CacheSize),1);

N_min_Zipf_RAND=zeros(length(CacheSize),1);
N_max_Zipf_RAND=zeros(length(CacheSize),1);
%% 3 class Uniform Distribution RAND (Random)

for nn=1:length(Prob_a)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    for cache=1:length(CacheSize)
        Probability_producers(1,:)=ProducersProbability_Uni(:,cache);
        N_a=CacheSize(cache);
        N_b=CacheSize(cache);
        N_c=Producers-N_a-N_b;
        Pop_producers=[N_a N_b N_c];
        Freshness_requirment=Freshness_Uni(:,cache);
        
        for jj=1:length(ProbForSavingVectorR1)
            ProbForSavingR1=ProbForSavingVectorR1(jj);
% memoryR1_LeastExpe and memoryR2_LeastExpe have following structure.
% First Column: Producer number
% Second Column: t_inst at which it is being fetched from producer

            memoryR1_RAND=zeros(CacheSize(cache),2);

            Router1_hit_count=zeros(Producers,1);

            count1=0;
            count2=0;

            message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f Probability_a=%f'...
                            ,CacheSize(cache),ProbForSavingR1,Prob_a(nn));
            h=msgbox(message);
            clear message

            N_min_temp=0;
            N_max_temp=0;
            display('3 class Uniform Distribution');
            for ii=1:length(time)
%                 display(t_inst);
                produ=producersRequest_Uni(ii,cache);
                t_inst=time(ii);
%                 Frshness=Freshness(t_inst,1);
                [N_min_temp,N_max_temp]=router1_RAND_plain_3class(produ,t_inst,...
                                                            ProbForSavingR1,N_min_temp,N_max_temp);

            end
            N_min(jj,cache)=N_min_temp;
            N_max(jj,cache)=N_max_temp;

            delete(h);
            clear('h');
        end        

        R1_hit_count_Uni_RAND(:,cache)=Router1_hit_count;
    end
    N_min_3class_RAND(:,nn)=N_min;
    N_max_3class_RAND(:,nn)=N_max;
end
toc

% delete(h1);
% clear('h1');
display('Done Uniform');

% requests_Uni(1:3,1);
% R1_hit_count_Uni(1:3,1)
% temp1(:,1)=R1_hit_count(1:3,1);
% temp2(:,1)=requests(1:3,1);
% display('Hit rate Simulation for  2 class Uniform Distribution')
% temp1./temp2
% temp2=repmat(requests_Uni,1,length(CacheSize));
hit_rate_Simul_Uni_RAND=R1_hit_count_Uni_RAND./requests_Uni;
clear temp1 temp2
hit_rate_total_Sim_Uni_RAND=sum(R1_hit_count_Uni_RAND)/NumberOfRequests;


%% LRU Variables ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
R1_hit_count_Uni_LRU=zeros(Producers,length(CacheSize));

R1_hit_count_Zipf_LRU=zeros(Producers,length(CacheSize));

N_min_3class_LRU=zeros(length(CacheSize),1);
N_max_3class_LRU=zeros(length(CacheSize),1);

N_min_Zipf_LRU=zeros(length(CacheSize),1);
N_max_Zipf_LRU=zeros(length(CacheSize),1);

%% 3 class Uniform Distribution LRU (Least Recently Used)
for nn=1:length(Prob_a)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    for cache=1:length(CacheSize)
        Probability_producers(1,:)=ProducersProbability_Uni(:,cache);
        N_a=CacheSize(cache);
        N_b=CacheSize(cache);
        N_c=Producers-N_a-N_b;
        Pop_producers=[N_a N_b N_c];
        Freshness_requirment=Freshness_Uni(:,cache);
        for jj=1:length(ProbForSavingVectorR1)
            ProbForSavingR1=ProbForSavingVectorR1(jj);
% memoryR1_LRU have following structure.
% First Column: t_stamp
% Second Column: Producer number
% Third Column: t_inst at which it is being fetched from producer
            memoryR1_LRU=zeros(CacheSize(cache),3);

            Router1_hit_count=zeros(Producers,1);

            message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f and Probability_a=%f'...
                            ,CacheSize(cache),ProbForSavingR1,Prob_a(nn));
            h=msgbox(message);
            clear message

            N_min_temp=0;
            N_max_temp=0;
            display('3 class Uniform Distribution');
            for ii=1:length(time)
%                 display(t_inst);
                produ=producersRequest_Uni(ii,cache);
                t_inst=time(ii);
%                 Frshness=Freshness(t_inst,1);
                [N_min_temp,N_max_temp]=router1_LRU_plain_3class(produ,t_inst,...
                                                            ProbForSavingR1,N_min_temp,N_max_temp);%,Sele_1,Sele_2);

            end
            N_min(jj,cache)=N_min_temp;
            N_max(jj,cache)=N_max_temp;

            delete(h);
            clear('h');
        end

        R1_hit_count_Uni_LRU(:,cache)=Router1_hit_count;
    end
    
    N_min_3class_LRU(:,nn)=N_min;
    N_max_3class_LRU(:,nn)=N_max;
end
toc


display('Done Uniform');

hit_rate_Simul_Uni_LRU=R1_hit_count_Uni_LRU./requests_Uni;
clear temp1 temp2
hit_rate_total_Sim_Uni_LRU=sum(R1_hit_count_Uni_LRU)/NumberOfRequests;

%% ###################################### Zipf Distribution with parameter beta #######################################
nn=1:Producers;
ProducersProbability_Zipf(1,:)=(nn.^-beta)/sum((nn.^-beta));

producersRequest_Zipf(:,1)=datasample(1:Producers,NumberOfRequests,'Weights',ProducersProbability_Zipf');
% temp1=repmat(ProducersProbability_Zipf',1,length(CacheSize));
% Freshness_Zipf=Freshness_const./temp1;
Freshness_Zipf=zeros(Producers,length(CacheSize));
for ii=1:length(CacheSize)
    Freshness_Zipf(:,ii)=sqrt((const*ProducersProbability_Zipf(1,1))./ProducersProbability_Zipf');
end

clear temp1;


%% Zipf distribution with parameter beta Least Expected
for nn=1:length(beta)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    Probability_producers(1,:)=ProducersProbability_Zipf(nn,:);

    for cache=1:length(CacheSize)
        N_a=CacheSize(cache);
        N_b=CacheSize(cache);
        N_c=Producers-N_a-N_b;
        Pop_producers=[N_a N_b N_c];
        Freshness_requirment=Freshness_Zipf(:,cache);

        for jj=1:length(ProbForSavingVectorR1)
            ProbForSavingR1=ProbForSavingVectorR1(jj);
    % memoryR1_LeastExpe and memoryR2_LeastExpe have following structure.
    % First Column: latest time_instant when data was being used under
    % condition it was fresh.
    % Second Column: Producer number
    % Third Column: time_stamp at ehich data for corresponding producer was
    % being fetched and stored.
            memoryR1_LeastExpe=zeros(CacheSize(cache),2);

            Router1_hit_count=zeros(Producers,1);

            count1=0;
            count2=0;

            message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f and beta=%f'...
                            ,CacheSize(cache),ProbForSavingR1,beta(nn));
            h=msgbox(message);
            clear message

            N_min_temp=0;
            N_max_temp=0;
            display('Zipf Distribution');
            for ii=1:length(time)
    %                 display(t_inst);
                produ=producersRequest_Zipf(ii,nn);
                t_inst=time(ii);
    %                 Frshness=Freshness(t_inst,1);
                [N_min_temp,N_max_temp]=router1_LeastExpe_plain_3class(produ,t_inst,...
                                                            ProbForSavingR1,N_min_temp,N_max_temp);

            end
            N_min(jj,cache)=N_min_temp; % jj-> row number; kk-> column number
            N_max(jj,cache)=N_max_temp;

            delete(h);
            clear('h');
        end

        R1_hit_count_Zipf_LeastExpe(:,cache)=Router1_hit_count;
    end
    N_min_Zipf_LeastExpe(:,nn)=N_min;
    N_max_Zipf_LeastExpe(:,nn)=N_max;
end
toc
display('Done Zipf');

requests_Zipf=zeros(Producers,length(CacheSize));
for ii=1:NumberOfRequests
    requests_Zipf(producersRequest_Zipf(ii,1),:)=requests_Zipf(producersRequest_Zipf(ii,1),:)+1;
end


% requests_Uni(1:3,1)
% R1_hit_count_Uni(1:3,1)
% temp1(:,1)=R1_hit_count(1:3,1);
% temp2(:,1)=requests(1:3,1);
% display('Hit rate Simulation for Zipf Distribution')
% temp1./temp2
% R1_hit_count_Zipf./requests_Zipf
% temp2=repmat(requests_Zipf,1,1,length(CacheSize));

hit_rate_Simul_Zipf_LeastExpe=R1_hit_count_Zipf_LeastExpe./requests_Zipf;
clear temp1 temp2
hit_rate_total_Sim_Zipf_LeastExpe=sum(R1_hit_count_Zipf_LeastExpe)/NumberOfRequests;

%% Zipf distribution with parameter beta SMP (Least Frequently Used)
for nn=1:length(beta)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    Probability_producers(1,:)=ProducersProbability_Zipf(nn,:);

    for cache=1:length(CacheSize)
        N_a=CacheSize(cache);
        N_b=CacheSize(cache);
        N_c=Producers-N_a-N_b;
        Pop_producers=[N_a N_b N_c];
        Freshness_requirment=Freshness_Zipf(:,cache);
        
        for jj=1:length(ProbForSavingVectorR1)
            ProbForSavingR1=ProbForSavingVectorR1(jj);
    % memoryR1_LeastExpe and memoryR2_LeastExpe have following structure.
    % First Column: latest time_instant when data was being used under
    % condition it was fresh.
    % Second Column: Producer number
    % Third Column: time_stamp at ehich data for corresponding producer was
    % being fetched and stored.
            memoryR1_SMP=zeros(CacheSize(cache),2);

            Router1_hit_count=zeros(Producers,1);

            count1=0;
            count2=0;

            message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f and beta=%f'...
                            ,CacheSize(cache),ProbForSavingR1,beta(nn));
            h=msgbox(message);
            clear message

            N_min_temp=0;
            N_max_temp=0;
            display('Zipf Distribution');
            for ii=1:length(time)
    %                 display(t_inst);
                produ=producersRequest_Zipf(ii,nn);
                t_inst=time(ii);
    %                 Frshness=Freshness(t_inst,1);
                [N_min_temp,N_max_temp]=router1_SMP_3class(produ,t_inst,...
                                                            ProbForSavingR1,N_min_temp,N_max_temp);

            end
            N_min(jj,cache)=N_min_temp; % jj-> row number; kk-> column number
            N_max(jj,cache)=N_max_temp;

            delete(h);
            clear('h');
        end

        R1_hit_count_Zipf_SMP(:,cache)=Router1_hit_count;
    end
    N_min_Zipf_SMP(:,nn)=N_min;
    N_max_Zipf_SMP(:,nn)=N_max;
end
toc
display('Done Zipf');


% requests_Uni(1:3,1)
% R1_hit_count_Uni(1:3,1)
% temp1(:,1)=R1_hit_count(1:3,1);
% temp2(:,1)=requests(1:3,1);
% display('Hit rate Simulation for Zipf Distribution')
% temp1./temp2
% R1_hit_count_Zipf./requests_Zipf
% temp2=repmat(requests_Zipf,1,1,length(CacheSize));

hit_rate_Simul_Zipf_SMP=R1_hit_count_Zipf_SMP./requests_Zipf;
clear temp1 temp2
hit_rate_total_Sim_Zipf_SMP=sum(R1_hit_count_Zipf_SMP)/NumberOfRequests;

%% Zipf distribution with parameter beta RAND (RANDOM)
for nn=1:length(beta)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    Probability_producers(1,:)=ProducersProbability_Zipf(nn,:);

    for cache=1:length(CacheSize)
        N_a=CacheSize(cache);
        N_b=CacheSize(cache);
        N_c=Producers-N_a-N_b;
        Pop_producers=[N_a N_b N_c];
        Freshness_requirment=Freshness_Zipf(:,cache);
        
        for jj=1:length(ProbForSavingVectorR1)
            ProbForSavingR1=ProbForSavingVectorR1(jj);
    % memoryR1_LeastExpe and memoryR2_LeastExpe have following structure.
    % First Column: latest time_instant when data was being used under
    % condition it was fresh.
    % Second Column: Producer number
    % Third Column: time_stamp at ehich data for corresponding producer was
    % being fetched and stored.
            memoryR1_RAND=zeros(CacheSize(cache),2);

            Router1_hit_count=zeros(Producers,1);

            count1=0;
            count2=0;

            message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f and beta=%f'...
                            ,CacheSize(cache),ProbForSavingR1,beta(nn));
            h=msgbox(message);
            clear message

            N_min_temp=0;
            N_max_temp=0;
            display('Zipf Distribution');
            for ii=1:length(time)
    %                 display(t_inst);
                produ=producersRequest_Zipf(ii,nn);
                t_inst=time(ii);
    %                 Frshness=Freshness(t_inst,1);
                [N_min_temp,N_max_temp]=router1_RAND_plain_3class(produ,t_inst,...
                                                            ProbForSavingR1,N_min_temp,N_max_temp);

            end
            N_min(jj,cache)=N_min_temp; % jj-> row number; kk-> column number
            N_max(jj,cache)=N_max_temp;

            delete(h);
            clear('h');
        end

        R1_hit_count_Zipf_RAND(:,cache)=Router1_hit_count;
    end
    N_min_Zipf_RAND(:,nn)=N_min;
    N_max_Zipf_RAND(:,nn)=N_max;
end
toc
display('Done Zipf');


% requests_Uni(1:3,1)
% R1_hit_count_Uni(1:3,1)
% temp1(:,1)=R1_hit_count(1:3,1);
% temp2(:,1)=requests(1:3,1);
% display('Hit rate Simulation for Zipf Distribution')
% temp1./temp2
% R1_hit_count_Zipf./requests_Zipf
% temp2=repmat(requests_Zipf,1,1,length(CacheSize));

hit_rate_Simul_Zipf_RAND=R1_hit_count_Zipf_RAND./requests_Zipf;
clear temp1 temp2
hit_rate_total_Sim_Zipf_RAND=sum(R1_hit_count_Zipf_RAND)/NumberOfRequests;

%% Zipf distribution with parameter beta LRU (Least Recently Used)
for nn=1:length(beta)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    Probability_producers(1,:)=ProducersProbability_Zipf(nn,:);

    for cache=1:length(CacheSize)
        N_a=CacheSize(cache);
        N_b=CacheSize(cache);
        N_c=Producers-N_a-N_b;
        Pop_producers=[N_a N_b N_c];
        Freshness_requirment=Freshness_Zipf(:,cache);
        
        for jj=1:length(ProbForSavingVectorR1)
            ProbForSavingR1=ProbForSavingVectorR1(jj);
% memoryR1_LRU and memoryR2_LRU have following structure.
% First Column: latest time_instant when data was being used under
% condition it was fresh.
% Second Column: Producer number
% Third Column: time_stamp at which data for corresponding producer was
% being fetched and stored.
            memoryR1_LRU=zeros(CacheSize(cache),3);

            Router1_hit_count=zeros(Producers,1);

            message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f and beta=%f'...
                            ,CacheSize(cache),ProbForSavingR1,beta(nn));
            h=msgbox(message);
            clear message

            N_min_temp=0;
            N_max_temp=0;
            display('Zipf Distribution');
            for ii=1:length(time)
%                 display(t_inst);
                produ=producersRequest_Zipf(ii,nn);
                t_inst=time(ii);
%                 Frshness=Freshness(t_inst,1);
                [N_min_temp,N_max_temp]=router1_LRU_plain_3class(produ,t_inst,...
                                                            ProbForSavingR1,N_min_temp,N_max_temp);%,Sele_1,Sele_2);

            end
            N_min(jj,cache)=N_min_temp;
            N_max(jj,cache)=N_max_temp;

            delete(h);
            clear('h');
        end
        
        R1_hit_count_Zipf_LRU(:,cache)=Router1_hit_count;
    end
    N_min_Zipf_LRU(:,nn)=N_min;
    N_max_Zipf_LRU(:,nn)=N_max;
end
toc
display('Done Zipf');

% temp2=repmat(requests_Zipf,1,1,length(CacheSize));
hit_rate_Simul_Zipf_LRU=R1_hit_count_Zipf_LRU./requests_Zipf;
clear temp1 temp2
hit_rate_total_Sim_Zipf_LRU=sum(R1_hit_count_Zipf_LRU)/NumberOfRequests;
%% Therotical Upper Bound 3 class Distribution
% Freshness=zeros(Producers,length(CacheSize));
% for cache=1:length(CacheSize)
%     N_a=CacheSize(cache);
%     N_b=CacheSize(cache);
%     N_c=Producers-2*CacheSize(cache);
%     Freshness(:,cache)=[repmat(F_a,N_a,1);repmat(F_b,N_b,1);repmat(F_c,N_c,1)];
% end

upperBound1_Uni=sum(((ProducersProbability_Uni.^2).*Freshness_Uni)./(ones(Producers,length(CacheSize))+ProducersProbability_Uni.*Freshness_Uni));
upperBound2_Uni=zeros(1,length(CacheSize));
for cache=1:length(CacheSize)
    upperBound2_Uni(1,cache)=sum(ProducersProbability_Uni(1:CacheSize(cache),cache));
end

upperBoundMin_Uni=min(upperBound1_Uni,upperBound2_Uni);

%% Therotical Upper Bound Zipf Distribution
temp1=repmat(ProducersProbability_Zipf',1,length(CacheSize));
upperBound1_Zipf=sum(((temp1.^2).*Freshness_Zipf)./(ones(Producers,length(CacheSize))+temp1.*Freshness_Zipf));

upperBound2_Zipf=zeros(1,length(CacheSize));
for cache=1:length(CacheSize)
    upperBound2_Zipf(1,cache)=sum(temp1(1:CacheSize(cache),cache));
end

upperBoundMin_Zipf=min(upperBound1_Zipf,upperBound2_Zipf);

%% Result Plot
% myplot(xinput,yinputMatrix,xlabel1,ylabel1,title1,legend1,saveFigAs)
clear temp1;
temp1=cd;
xinput(:,1)=CacheSize;
yinputMatrix=horzcat(upperBoundMin_Uni',hit_rate_total_Sim_Uni_LeastExpe',hit_rate_total_Sim_Uni_LRU',hit_rate_total_Sim_Uni_RAND',hit_rate_total_Sim_Uni_SMP');
xlabel1=sprintf('Cache size');
ylabel1=sprintf('Hit rate');
% title1=sprintf('Hit rate (p_{hit}) Vs Cache size');
directory='D:\IoT\IoT\31Jan\LeastExpected\New_Results24092017';
legend1={sprintf('UB'),sprintf('LU'),sprintf('LRU'),sprintf('RAND'),sprintf('SMP')};
saveFigAs=sprintf('Hit_rate_Vs_Cache_Size_policies_Uniform');
myplot(xinput,yinputMatrix,xlabel1,ylabel1,legend1,saveFigAs,directory);
cd(temp1);

yinputMatrix=horzcat(upperBoundMin_Zipf',hit_rate_total_Sim_Zipf_LeastExpe',hit_rate_total_Sim_Zipf_LRU',hit_rate_total_Sim_Zipf_RAND',hit_rate_total_Sim_Zipf_SMP');
xlabel1=sprintf('Cache size');
ylabel1=sprintf('Hit rate');
% title1=sprintf('Hit rate (p_{hit}) Vs Cache size');
directory='D:\IoT\IoT\31Jan\LeastExpected\New_Results24092017';
legend1={sprintf('UB'),sprintf('LU'),sprintf('LRU'),sprintf('RAND'),sprintf('SMP')};
saveFigAs=sprintf('Hit_rate_Vs_Cache_Size_policies_Zipf');
myplot(xinput,yinputMatrix,xlabel1,ylabel1,legend1,saveFigAs,directory);
%% always change the dataname for saving. Keep it simple and discriptive.
% temp1=cd;
cd('D:\IoT\IoT\31Jan\LeastExpected\New_Results24092017\Data')
save('cmp_cache_policies_3Class_Sqrt_Freshness');
