function results = Simulator
  
  % settings
  N = 100;                             % pop size
  graphType = 'Erdos-Renyi';           % initial network structure
  numRounds = 2000;                    % total number of rounds
  payoffMatrix = [3, 0; 5, 1];         % payoff matrix for the game
  pCooperator = 0.5;                   % initial probability of playing allC
  pRewireRound = 0.1;                  % probability of a rewire round
  pMutation = 0.1;                     % probability of mutating on play round
  luceExponentSD = 2;                  % standard deviation of Luce choice exp
  
  % derivative parameters
  isRewireRound = rand(1,numRounds) < pRewireRound;

  % create the initial pop
  pop.graph = MakeAdjacencyMatrix(graphType, N);
  pop.strategies = generateRandomStrategies(N); %
  
  for i = 1:numRounds
    
    % rewire round
    if(rand < pRewireRound)
      %display('rewire')
      % pick a random connection, player 1 makes the choice
      players = shuffle(1:N);
      player1strategy = pop.strategies(players(1),:);
      player2strategy = pop.strategies(players(2),:);
      rewire = rand() < player1strategy(player2strategy(1)+1);
      if(rewire)
        degree = full(sum(pop.graph > 0, 2));
        unlinkedIndices = find(~pop.graph(players(1),:));    
        L = player1strategy(4);
        p = (degree(unlinkedIndices).^L)/sum(degree(unlinkedIndices).^L);
        newConnection = unlinkedIndices(randp(p));
        
        % remove the old connection; form the new one
        pop.graph(players(1),players(2)) = 0;
        pop.graph(players(1),newConnection) = 1;
        pop.graph = rownormalize(pop.graph);
      end

    % playing round  
    else
      %display('--play--')
      if(rand() < pMutation)
        pop.strategies(randi(N),:) = generateRandomStrategies(1);
      else
        % simulate a pairwise comparison process: choose a pair of players,
        % assign one the role of student, the other the role of teacher.
        % the student takes on the teacher's strategy with probability
        % proportional to payoff. implicit in this model is a replacement
        % network that is complete.
        players = shuffle(1:N);
        payoffs(1) = playWithNeighbors(players(1)); % student
        payoffs(2) = playWithNeighbors(players(2)); % teacher
        
        isReplace = rand() < payoffs(2)/sum(payoffs);
        if(isReplace)
          pop.strategies(players(1),:) = pop.strategies(players(2),:);
        end
      end
    end
  end
  results = -1;
  
  %
  function payoff = playWithNeighbors(player)
    neighbors = find(pop.graph(player,:));
    for k = 1:length(neighbors)
      payoff(k) = playGame(pop.strategies(player,:), ...
                           pop.strategies(neighbors(k),:));
    end
    payoff = sum(payoff);
  end
  
  % plays a game between players with strategies s1 and s2
  function [p1, p2] = playGame(s1, s2)
    p1 = payoffMatrix(s1(1),s2(1));
    p2 = payoffMatrix(s2(1),s1(1));
  end
  
  % TODO: add noise epsilon
  function s = generateRandomStrategies(N)
    s(:,1) = 1 + (rand(N,1) > pCooperator); % allC = 1, allD = 2
    s(:,2) = rand(N,1);                     % P(rewire|opponentIsCooperator)
    s(:,3) = rand(N,1);                     % P(rewire|opponentIsDefector)
    s(:,4) = luceExponentSD*randn(N,1);     % Luce choice exponent
  end
end
