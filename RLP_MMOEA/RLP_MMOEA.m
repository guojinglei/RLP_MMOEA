classdef RLP_MMOEA < ALGORITHM
% <multi> <real/integer> <multimodal>
% K --- 4 --- Parameter for LPA
% O --- 5 --- Parameter for Insert
% theta --- 0.2 --- Parameter for Outer Point
% garma --- 1.2 --- Parameter for Outer Point


    methods
        function main(Algorithm,Problem)

            [K, O, theta, garma] = Algorithm.ParameterSet(4, 5, 0.2, 1.2);

            G = ceil(Problem.maxFE/(2 * Problem.N));
            gen = 1;
            %% Generate random population
            Population1 = Problem.Initialization();
            Population2 = Problem.Initialization();
            [Fitness1,D_Dec,D_Pop] = CalFitness(Population1.objs,Population1.decs);
            %% Optimization
            while Algorithm.NotTerminated(Population1)
                
                EvoState = gen / G;
                CurEps = -(1/(1-pi.^EvoState))-0.5; % Compute the current epsilon value
  
                MatingPool1 = TournamentSelection(2,Problem.N,D_Dec,D_Pop,Fitness1);
                Offspring1  = OperatorGAhalf(Problem,Population1(MatingPool1));
                
                %% LPA Clustering Method
                [~,adjacent_matrix] = sort(pdist2(Population2.decs,Population2.decs),2);
                adjacent_matrix = adjacent_matrix(:,2:K+1); 
                label = 1:1:Problem.N;
               [ Labelnew ] = LPA( adjacent_matrix,label );
                
                %% Generate Offspring
                Population_new = [];
                
                ulabel = unique(Labelnew);
                C = length(ulabel);
                
                for i = 1:C
                    cluster = find(Labelnew==ulabel(i));
                    [Fitness2,D] = CalFitnessDecEpsilon(Population2(cluster).objs,Population2(cluster).decs,CurEps);
                    MatingPool2 = TournamentSelection(2,2*length(cluster),D,Fitness2);
                    MatingPop = Population2(cluster);
                    Offspring2  = OperatorDE(Problem,MatingPop,MatingPop(MatingPool2(1:end/2)),MatingPop(MatingPool2(end/2+1:end)));
                    Population_new = [Population_new,Offspring2];
                end
                
                Offspring2 = Population_new;
                %% Environmental Selection
                [Population1,Fitness1,D_Dec,D_Pop] = EnvironmentalSelection([Population1,Offspring1,Offspring2],Problem.N);
                [Population2,~,~] = EnvironmentalSelectionDec([Population2,Offspring1,Offspring2],Problem.N,CurEps);
                    
                gen = ceil(Problem.FE/(2 * Problem.N));


                if EvoState > theta 
                    dist = pdist2(Population1.decs, Population1.decs); 
                    [dist, ~] = sort(dist, 2);
                    for i = 1:Problem.N
                        d(i,:) = dist(i, 2 : K + 1);
                        outdist(i) =  sum(d(i,:)) / K;
                    end
                    avg_outdist = sum(outdist) / Problem.N;

                    Pop_temp = [];
                    for i = 1:Problem.N
                        if (outdist(i) > garma * avg_outdist)
                            Pop_temp = [Pop_temp, Population1(i)];
                        end
                    end

                    if ~isempty(Pop_temp)
                        ddist = pdist2(Pop_temp.decs, Population1.decs);
        
                        [~, index] = sort(ddist, 2);
                        len = size(index, 1);
                        Off = [];
                        for i = 1 : len 
                            Off = [Off; Insert(Population1(index(i, 1:O+1)))];
                        end
                        [Off,~,~] = unique(Off, "rows");
                        Off_new = Problem.Evaluation(Off);
        
                        [Population1,Fitness1,D_Dec,D_Pop] = EnvironmentalSelection([Population1,Off_new],Problem.N);
                        [Population2,~,~] = EnvironmentalSelectionDec([Population2,Off_new],Problem.N,CurEps);
                    end           
                end
                

            end
        end
    end
end