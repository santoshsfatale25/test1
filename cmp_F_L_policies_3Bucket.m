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
% close all;
% clc;
%%
NumberOfIterations=10^2;
NumberOfRequests=10^5;
Producers=100; % Number of Producers
global Pop_producers

global Freshness_requirment
const=100;
F_a=[5 10 20 50 100];
F_b=const*F_a;
F_c=F_a;
% Freshness_requirment=[F_a F_b F_c];

global Router1_hit_count

ProbForSavingVectorR1=1;%0.2:0.2:1.0;%1.0;
CacheSize=20;%10:5:40;

Prob_a=0.4;%0.25:0.05:0.45;
beta=0.8;%0.5:0.3:1.7;
hit_rate_total_Sim_Uni_LeastExpe=zeros(NumberOfIterations,length(F_a));
hit_rate_total_Sim_Uni_SMP=zeros(NumberOfIterations,length(F_a));
hit_rate_total_Sim_Uni_RAND=zeros(NumberOfIterations,length(F_a));
hit_rate_total_Sim_Uni_LRU=zeros(NumberOfIterations,length(F_a));
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

% Characteristic_time_Uni=zeros(Producers,length(CacheSize));
% Characteristic_time_Zipf=zeros(Producers,length(CacheSize));

global memoryR1_LeastExpe memoryR1_LRU memoryR1_SMP memoryR1_RAND Probability_producers

global count1 count2 % Checks cache is empty or not.

%% ######################################### 3 class Uniform Distribution ################################################

tic;
Prob_b=0.99*Prob_a;
Prob_c=ones(1,length(Prob_a))-Prob_a-Prob_b;

N_a=CacheSize;
N_b=CacheSize;
N_c=Producers-N_a-N_b;
ProducersProbability_Uni(:,1)=[repmat(Prob_a./N_a,N_a,1);repmat(Prob_b./N_b,N_b,1);repmat(Prob_c./N_c,N_c,1)];
Freshness_Uni=[repmat(F_a,N_a,1);repmat(F_b,N_b,1);repmat(F_c,N_c,1)];

for kk=1:NumberOfIterations
    display(sprintf('Iteration Number=%d',kk));
    %% Exponential inter-arrival time
    time=cumsum(exprnd(1,NumberOfRequests,1));
%    producersRequest_Uni=zeros(Producers,1);
%    for pp=1:length(Producers)
%       r=rand(1,NumberOfRequests);
%       for ii=1:NumberOfRequests
%           temp1(1,:)=ProducersProbability_Uni;
%           producersRequest_Uni(ii,1) = sum(r(1,ii) >= cumsum([0, temp1]));
%       end
%    end
    clear r temp1;
   producersRequest_Uni(:,1)=datasample(1:Producers,NumberOfRequests,'Weights',ProducersProbability_Uni);
    % temp1=repmat(ProducersProbability_Uni,1,length(F_a));
    % Freshness_const=Freshness_Uni.*temp1;

    Probability_producers(1,:)=ProducersProbability_Uni(:,1);

    %% Least Expected Variables :::::::::::::::::::::::::::::::::::::::::::::::
    R1_hit_count_Uni_LeastExpe=zeros(Producers,length(F_a));

    R1_hit_count_Zipf_LeastExpe=zeros(Producers,length(F_a));

    N_min_3class_LeastExpe=zeros(length(CacheSize),length(F_a));
    N_max_3class_LeastExpe=zeros(length(CacheSize),length(F_a));

    N_min_Zipf_LeastExpe=zeros(length(CacheSize),length(F_a));
    N_max_Zipf_LeastExpe=zeros(length(CacheSize),length(F_a));

    %% 3 class Uniform Distribution Least Expected

    for nn=1:length(F_a)
        N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
        N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
        for cache=1:length(CacheSize)
    %         Probability_producers(1,:)=ProducersProbability_Uni(:,nn);
            N_a=CacheSize(cache);
            N_b=CacheSize(cache);
            N_c=Producers-N_a-N_b;
            Pop_producers=[N_a N_b N_c];
            Freshness_requirment=Freshness_Uni(:,nn);

            for jj=1:length(ProbForSavingVectorR1)
                ProbForSavingR1=ProbForSavingVectorR1(jj);
    % memoryR1_LeastExpe and memoryR2_LeastExpe have following structure.
    % First Column: Producer number
    % Second Column: t_inst at which it is being fetched from producer

                memoryR1_LeastExpe=zeros(CacheSize(cache),2);

                Router1_hit_count=zeros(Producers,1);

                count1=0;
                count2=0;

                message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f Probability_a=%f'...
                                ,CacheSize(cache),ProbForSavingR1,Prob_a);
                h=msgbox(message);
                clear message

                N_min_temp=0;
                N_max_temp=0;
                display('3 class Uniform Distribution');
                for ii=1:length(time)
    %                 display(t_inst);
                    produ=producersRequest_Uni(ii,1);
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

        end
        R1_hit_count_Uni_LeastExpe(:,nn)=Router1_hit_count;
        N_min_3class_LeastExpe(:,nn)=N_min;
        N_max_3class_LeastExpe(:,nn)=N_max;

    end
    toc

    % delete(h1);
    % clear('h1');
    display('Done Uniform');


    requests_Uni=zeros(Producers,1);
    for ii=1:NumberOfRequests
        requests_Uni(producersRequest_Uni(ii,1),1)=requests_Uni(producersRequest_Uni(ii,1),1)+1;
    end



    % requests_Uni(1:3,1);
    % R1_hit_count_Uni(1:3,1)
    % temp1(:,1)=R1_hit_count(1:3,1);
    % temp2(:,1)=requests(1:3,1);
    % display('Hit rate Simulation for  2 class Uniform Distribution')
    % temp1./temp2
    % temp2=repmat(requests_Uni,1,length(F_a));
    % hit_rate_Simul_Uni_LeastExpe=R1_hit_count_Uni_LeastExpe./temp2;
    % clear temp1 temp2
    hit_rate_total_Sim_Uni_LeastExpe(kk,:)=sum(R1_hit_count_Uni_LeastExpe)/NumberOfRequests;


    %% SMP Variables ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    R1_hit_count_Uni_SMP=zeros(Producers,length(F_a));

    R1_hit_count_Zipf_SMP=zeros(Producers,length(F_a));

    N_min_3class_SMP=zeros(length(CacheSize),length(F_a));
    N_max_3class_SMP=zeros(length(CacheSize),length(F_a));

    N_min_Zipf_SMP=zeros(length(CacheSize),length(F_a));
    N_max_Zipf_SMP=zeros(length(CacheSize),length(F_a));

    %% 3 class Uniform Distribution SMP (Least Frequently Used)

    for nn=1:length(F_a)
        N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
        N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
        for cache=1:length(CacheSize)
    %         Probability_producers(1,:)=ProducersProbability_Uni(:,nn);
            N_a=CacheSize(cache);
            N_b=CacheSize(cache);
            N_c=Producers-N_a-N_b;
            Pop_producers=[N_a N_b N_c];
            Freshness_requirment=Freshness_Uni(:,nn);

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
                                ,CacheSize(cache),ProbForSavingR1,Prob_a);
                h=msgbox(message);
                clear message

                N_min_temp=0;
                N_max_temp=0;
                display('3 class Uniform Distribution');
                for ii=1:length(time)
    %                 display(t_inst);
                    produ=producersRequest_Uni(ii,1);
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

        end
        R1_hit_count_Uni_SMP(:,nn)=Router1_hit_count;
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
    % temp2=repmat(requests_Uni,1,length(F_a));
    % hit_rate_Simul_Uni_SMP=R1_hit_count_Uni_SMP./temp2;
    % clear temp1 temp2
    hit_rate_total_Sim_Uni_SMP(kk,:)=sum(R1_hit_count_Uni_SMP)/NumberOfRequests;


    %% RAND Variables ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    R1_hit_count_Uni_RAND=zeros(Producers,length(F_a));

    R1_hit_count_Zipf_RAND=zeros(Producers,length(F_a));

    N_min_3class_RAND=zeros(length(CacheSize),length(F_a));
    N_max_3class_RAND=zeros(length(CacheSize),length(F_a));

    N_min_Zipf_RAND=zeros(length(CacheSize),length(F_a));
    N_max_Zipf_RAND=zeros(length(CacheSize),length(F_a));
    %% 3 class Uniform Distribution RAND (Random)

    for nn=1:length(F_a)
        N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
        N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
        for cache=1:length(CacheSize)
    %         Probability_producers(1,:)=ProducersProbability_Uni(:,nn);
            N_a=CacheSize(cache);
            N_b=CacheSize(cache);
            N_c=Producers-N_a-N_b;
            Pop_producers=[N_a N_b N_c];
            Freshness_requirment=Freshness_Uni(:,nn);

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
                                ,CacheSize(cache),ProbForSavingR1,Prob_a);
                h=msgbox(message);
                clear message

                N_min_temp=0;
                N_max_temp=0;
                display('3 class Uniform Distribution');
                for ii=1:length(time)
    %                 display(t_inst);
                    produ=producersRequest_Uni(ii,1);
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

        end
        R1_hit_count_Uni_RAND(:,nn)=Router1_hit_count;
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
    % temp2=repmat(requests_Uni,1,length(F_a));
    % hit_rate_Simul_Uni_RAND=R1_hit_count_Uni_RAND./temp2;
    % clear temp1 temp2
    hit_rate_total_Sim_Uni_RAND(kk,:)=sum(R1_hit_count_Uni_RAND)/NumberOfRequests;


    %% LRU Variables ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    R1_hit_count_Uni_LRU=zeros(Producers,length(F_a));

    R1_hit_count_Zipf_LRU=zeros(Producers,length(F_a));

    N_min_3class_LRU=zeros(length(CacheSize),length(F_a));
    N_max_3class_LRU=zeros(length(CacheSize),length(F_a));

    N_min_Zipf_LRU=zeros(length(CacheSize),length(F_a));
    N_max_Zipf_LRU=zeros(length(CacheSize),length(F_a));

    %% 3 class Uniform Distribution LRU (Least Recently Used)

    for nn=1:length(F_a)
        N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
        N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
        for cache=1:length(CacheSize)
    %         Probability_producers(1,:)=ProducersProbability_Uni(:,nn);
            N_a=CacheSize(cache);
            N_b=CacheSize(cache);
            N_c=Producers-N_a-N_b;
            Pop_producers=[N_a N_b N_c];
            Freshness_requirment=Freshness_Uni(:,nn);

            for jj=1:length(ProbForSavingVectorR1)
                ProbForSavingR1=ProbForSavingVectorR1(jj);
    % memoryR1_LRU have following structure.
    % First Column: t_stamp
    % Second Column: Producer number
    % Third Column: t_inst at which it is being fetched from producer
                memoryR1_LRU=zeros(CacheSize(cache),3);

                Router1_hit_count=zeros(Producers,1);

                message=sprintf('Running for Cache Size=%d and ProbForSavingR1=%f and Probability_a=%f'...
                                ,CacheSize(cache),ProbForSavingR1,Prob_a);
                h=msgbox(message);
                clear message

                N_min_temp=0;
                N_max_temp=0;
                display('3 class Uniform Distribution');
                for ii=1:length(time)
    %                 display(t_inst);
                    produ=producersRequest_Uni(ii,1);
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

        end
        R1_hit_count_Uni_LRU(:,nn)=Router1_hit_count;
        N_min_3class_LRU(:,nn)=N_min;
        N_max_3class_LRU(:,nn)=N_max;
    end
    toc


    display('Done Uniform');

    % temp2=repmat(requests_Uni,1,length(F_a));
    % hit_rate_Simul_Uni_LRU=R1_hit_count_Uni_LRU./temp2;
    % clear temp1 temp2
    hit_rate_total_Sim_Uni_LRU(kk,:)=sum(R1_hit_count_Uni_LRU)/NumberOfRequests;
end
hit_rate_total_Sim_Uni_LeastExpe_Average=mean(hit_rate_total_Sim_Uni_LeastExpe);
hit_rate_total_Sim_Uni_LeastExpe_stdDev=std(hit_rate_total_Sim_Uni_LeastExpe,1,1);

hit_rate_total_Sim_Uni_LRU_Average=mean(hit_rate_total_Sim_Uni_LRU);
hit_rate_total_Sim_Uni_LRU_stdDev=std(hit_rate_total_Sim_Uni_LRU,1,1);

hit_rate_total_Sim_Uni_RAND_Average=mean(hit_rate_total_Sim_Uni_RAND);
hit_rate_total_Sim_Uni_RAND_stdDev=std(hit_rate_total_Sim_Uni_RAND,1,1);

hit_rate_total_Sim_Uni_SMP_Average=mean(hit_rate_total_Sim_Uni_SMP);
hit_rate_total_Sim_Uni_SMP_stdDev=std(hit_rate_total_Sim_Uni_SMP,1,1);

clear hit_rate_Simul_Uni_LeastExpe hit_rate_Simul_Uni_LRU hit_rate_Simul_Uni_RAND hit_rate_Simul_Uni_SMP
clear R1_hit_count_Uni_LeastExpe R1_hit_count_Uni_LRU R1_hit_count_Uni_RAND R1_hit_count_Uni_SMP
clear N_min_3class_LeastExpe N_min_3class_LRU N_min_3class_RAND N_min_3class_SMP N_max_3class_LeastExpe N_max_3class_LRU N_max_3class_RAND N_max_3class_SMP
clear produ t_inst N_min N_max N_min_temp N_max_temp

%% Therotical Upper Bound 3 class Distribution

temp1=repmat(ProducersProbability_Uni,1,length(F_a));
upperBound1_Uni=sum(((temp1.^2).*Freshness_Uni)./(ones(Producers,length(F_a))+temp1.*Freshness_Uni));

upperBound2_Uni=zeros(1,length(F_a));
for nn=1:length(F_a)
    upperBound2_Uni(1,nn)=sum(temp1(1:CacheSize,nn));
end
clear temp1

upperBoundMin_Uni=min(upperBound1_Uni,upperBound2_Uni);

%% Result Plot
% myplotNew(xinput,yinputMatrix_avg,yinputMatrix_stdDev,xlabel1,ylabel1,title1,legend1,xlim1,ylim1,saveFigAs,directory)
clear temp1;
temp1=cd;
xinput(:,1)=F_a;
yinputMatrix_avg=horzcat(upperBoundMin_Uni',hit_rate_total_Sim_Uni_LeastExpe_Average',hit_rate_total_Sim_Uni_LRU_Average',hit_rate_total_Sim_Uni_RAND_Average',hit_rate_total_Sim_Uni_SMP_Average');
yinputMatrix_stdDev=horzcat(hit_rate_total_Sim_Uni_LeastExpe_stdDev',hit_rate_total_Sim_Uni_LRU_stdDev',hit_rate_total_Sim_Uni_RAND_stdDev',hit_rate_total_Sim_Uni_SMP_stdDev');
xlabel1=sprintf('F_a');
ylabel1=sprintf('Cache hit ratio');
% title1=sprintf('Hit rate (p_{hit}) Vs Cache size');
directory='D:\IoT\IoT\31Jan\LeastExpected\NCC_Results';
legend1={sprintf('UB'),sprintf('LU'),sprintf('LRU'),sprintf('RAND'),sprintf('SMP')};
saveFigAs=sprintf('Hit_rate_Vs_F_a_policies_Uniform');
xlim1=[min(F_a) max(F_a)];
ylim1=[0 0.45];
myplotNew(xinput,yinputMatrix_avg,yinputMatrix_stdDev,xlabel1,ylabel1,legend1,xlim1,ylim1,saveFigAs,directory);
cd(temp1);

%% always change the dataname for saving. Keep it simple and discriptive.
% temp1=cd;
cd('D:\IoT\IoT\31Jan\LeastExpected\NCC_Results')
save('cmp_F_a_policies_3Class');





%% ###################################### Zipf Distribution with parameter beta #######################################
nn=1:Producers;
ProducersProbability_Zipf(1,:)=(nn.^-beta)/sum((nn.^-beta));
temp1=repmat(ProducersProbability_Zipf',1,length(F_a));
Freshness_Zipf=Freshness_const./temp1;

producersRequest_Zipf(:,1)=datasample(1:Producers,NumberOfRequests,'Weights',ProducersProbability_Zipf');

Probability_producers(1,:)=ProducersProbability_Zipf(1,:);

%% Zipf distribution with parameter beta Least Expected
for nn=1:length(F_a)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));

    for cache=1:length(CacheSize)
        N_a=CacheSize(cache);
        N_b=CacheSize(cache);
        N_c=Producers-N_a-N_b;
        Pop_producers=[N_a N_b N_c];
        Freshness_requirment=Freshness_Zipf(:,nn);

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
                            ,CacheSize(cache),ProbForSavingR1,beta);
            h=msgbox(message);
            clear message

            N_min_temp=0;
            N_max_temp=0;
            display('Zipf Distribution');
            for ii=1:length(time)
    %                 display(t_inst);
                produ=producersRequest_Zipf(ii,1);
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

    end
    R1_hit_count_Zipf_LeastExpe(:,nn)=Router1_hit_count;
    N_min_Zipf_LeastExpe(:,nn)=N_min;
    N_max_Zipf_LeastExpe(:,nn)=N_max;
end
toc
display('Done Zipf');

requests_Zipf=zeros(Producers,1);
for ii=1:NumberOfRequests
    requests_Zipf(producersRequest_Zipf(ii,1),1)=requests_Zipf(producersRequest_Zipf(ii,1),1)+1;
end


% requests_Uni(1:3,1)
% R1_hit_count_Uni(1:3,1)
% temp1(:,1)=R1_hit_count(1:3,1);
% temp2(:,1)=requests(1:3,1);
% display('Hit rate Simulation for Zipf Distribution')
% temp1./temp2
% R1_hit_count_Zipf./requests_Zipf
temp2=repmat(requests_Zipf,1,length(F_a));
hit_rate_Simul_Zipf_LeastExpe=R1_hit_count_Zipf_LeastExpe./temp2;
clear temp1 temp2
hit_rate_total_Sim_Zipf_LeastExpe=sum(R1_hit_count_Zipf_LeastExpe)/NumberOfRequests;

%% Zipf distribution with parameter beta SMP (Least Frequently Used)
for nn=1:length(F_a)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
%     Probability_producers(1,:)=ProducersProbability_Zipf(nn,:);

    for cache=1:length(CacheSize)
        N_a=CacheSize(cache);
        N_b=CacheSize(cache);
        N_c=Producers-N_a-N_b;
        Pop_producers=[N_a N_b N_c];
        Freshness_requirment=Freshness_Zipf(:,nn);

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
                            ,CacheSize(cache),ProbForSavingR1,beta);
            h=msgbox(message);
            clear message

            N_min_temp=0;
            N_max_temp=0;
            display('Zipf Distribution');
            for ii=1:length(time)
    %                 display(t_inst);
                produ=producersRequest_Zipf(ii,1);
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

    end
    R1_hit_count_Zipf_SMP(:,nn)=Router1_hit_count;    
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
temp2=repmat(requests_Zipf,1,length(F_a));
hit_rate_Simul_Zipf_SMP=R1_hit_count_Zipf_SMP./temp2;
clear temp1 temp2
hit_rate_total_Sim_Zipf_SMP=sum(R1_hit_count_Zipf_SMP)/NumberOfRequests;

%% Zipf distribution with parameter beta RAND (RANDOM)
for nn=1:length(F_a)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
%     Probability_producers(1,:)=ProducersProbability_Zipf(nn,:);

    for cache=1:length(CacheSize)
        N_a=CacheSize(cache);
        N_b=CacheSize(cache);
        N_c=Producers-N_a-N_b;
        Pop_producers=[N_a N_b N_c];
        Freshness_requirment=Freshness_Zipf(:,nn);
        
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
                            ,CacheSize(cache),ProbForSavingR1,beta);
            h=msgbox(message);
            clear message

            N_min_temp=0;
            N_max_temp=0;
            display('Zipf Distribution');
            for ii=1:length(time)
    %                 display(t_inst);
                produ=producersRequest_Zipf(ii,1);
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

    end
    R1_hit_count_Zipf_RAND(:,nn)=Router1_hit_count;
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
temp2=repmat(requests_Zipf,1,length(F_a));
hit_rate_Simul_Zipf_RAND=R1_hit_count_Zipf_RAND./temp2;
clear temp1 temp2
hit_rate_total_Sim_Zipf_RAND=sum(R1_hit_count_Zipf_RAND)/NumberOfRequests;

%% Zipf distribution with parameter beta LRU (Least Recently Used)
for nn=1:length(F_a)
    N_min=zeros(length(ProbForSavingVectorR1),length(CacheSize));
    N_max=zeros(length(ProbForSavingVectorR1),length(CacheSize));
%     Probability_producers(1,:)=ProducersProbability_Zipf(nn,:);

    for cache=1:length(CacheSize)
        N_a=CacheSize(cache);
        N_b=CacheSize(cache);
        N_c=Producers-N_a-N_b;
        Pop_producers=[N_a N_b N_c];
        Freshness_requirment=Freshness_Zipf(:,nn);
        
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
                            ,CacheSize(cache),ProbForSavingR1,beta);
            h=msgbox(message);
            clear message

            N_min_temp=0;
            N_max_temp=0;
            display('Zipf Distribution');
            for ii=1:length(time)
%                 display(t_inst);
                produ=producersRequest_Zipf(ii,1);
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
        
    end
    R1_hit_count_Zipf_LRU(:,nn)=Router1_hit_count;
    N_min_Zipf_LRU(:,nn)=N_min;
    N_max_Zipf_LRU(:,nn)=N_max;
end
toc
display('Done Zipf');

temp2=repmat(requests_Zipf,1,length(F_a));
hit_rate_Simul_Zipf_LRU=R1_hit_count_Zipf_LRU./temp2;
clear temp1 temp2
hit_rate_total_Sim_Zipf_LRU=sum(R1_hit_count_Zipf_LRU)/NumberOfRequests;


%% Therotical Upper Bound Zipf Distribution

temp1=repmat(ProducersProbability_Zipf',1,length(F_a));
upperBound1_Zipf=sum(((temp1.^2).*Freshness_Zipf)./(ones(Producers,length(F_a))+temp1.*Freshness_Zipf));

upperBound2_Zipf=zeros(1,length(F_a));
for nn=1:length(F_a)
    upperBound2_Zipf(1,nn)=sum(temp1(1:CacheSize,nn));
end
clear temp1

upperBoundMin_Zipf=min(upperBound1_Zipf,upperBound2_Zipf);

%% Result Plot
% myplot(xinput,yinputMatrix,xlabel1,ylabel1,title1,legend1,saveFigAs)
xinput(:,1)=F_a;
yinputMatrix=horzcat(upperBoundMin_Uni',hit_rate_total_Sim_Uni_LeastExpe',hit_rate_total_Sim_Uni_LRU',hit_rate_total_Sim_Uni_RAND',hit_rate_total_Sim_Uni_SMP');
xlabel1=sprintf('F_{low}');
ylabel1=sprintf('Hit rate');
% title1=sprintf('Hit rate (p_{hit}) Vs F_{low}');
legend1={sprintf('UB'),sprintf('LU'),sprintf('LRU'),sprintf('RAND'),sprintf('SMP')};
saveFigAs=sprintf('Hit_rate_Vs_F_L_policies_Uniform');
myplot(xinput,yinputMatrix,xlabel1,ylabel1,legend1,saveFigAs);
%% Result Plot

yinputMatrix=horzcat(upperBoundMin_Zipf',hit_rate_total_Sim_Zipf_LeastExpe',hit_rate_total_Sim_Zipf_LRU',hit_rate_total_Sim_Zipf_RAND',hit_rate_total_Sim_Zipf_SMP');
xlabel1=sprintf('F_{low}');
ylabel1=sprintf('Hit rate');
% title1=sprintf('Hit rate (p_{hit}) Vs F_{low}');
legend1={sprintf('UB'),sprintf('LU'),sprintf('LRU'),sprintf('RAND'),sprintf('SMP')};
saveFigAs=sprintf('Hit_rate_Vs_F_L_policies_Zipf');
myplot(xinput,yinputMatrix,xlabel1,ylabel1,legend1,saveFigAs);

%% always change the dataname for saving. Keep it simple and discriptive.
% temp1=cd;
cd('D:\IoT\IoT\31Jan\Least Expected\Data')
save('cmp_F_L_policies_3class'); % F_L is freshenss requirement for class 1 and class 3. 
