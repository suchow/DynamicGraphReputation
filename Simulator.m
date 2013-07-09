function results = Simulator(s)
  
  % derivative parameters
  isRewireRound = rand(1,s.numRounds) < s.pRewireRound;

  % create the initial pop
  pop.graph = MakeAdjacencyMatrix(s.graphType, s.N);
  pop.strategies = generateRandomStrategies(s.N); %
  
  for i = 1:s.numRounds
    
    % rewire round
    if(rand < s.pRewireRound)
      %display('rewire')
      % pick a random connection, player 1 makes the choice
      players = shuffle(1:s.N);
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
      if(rand() < s.pMutation)
        pop.strategies(randi(s.N),:) = generateRandomStrategies(1);
      else
        % simulate a pairwise comparison process: choose a pair of players,
        % assign one the role of student, the other the role of teacher.
        % the student takes on the teacher's strategy with probability
        % proportional to payoff. implicit in this model is a replacement
        % network that is complete.
        players = shuffle(1:s.N);
        payoffs(1) = playWithNeighbors(players(1)); % student
        payoffs(2) = playWithNeighbors(players(2)); % teacher
        
        isReplace = rand() < payoffs(2)/sum(payoffs);
        if(isReplace)
          pop.strategies(players(1),:) = pop.strategies(players(2),:);
        end
      end
    end
    % store results
    results.history(i) = pop;
  end
  
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
    p1 = s.payoffMatrix(s1(1),s2(1));
    p2 = s.payoffMatrix(s2(1),s1(1));
  end
  
  % TODO: add noise epsilon
  function str = generateRandomStrategies(N)
    str(:,1) = 1 + (rand(N,1) > s.pCooperator);% allC = 1, allD = 2
    str(:,2) = rand(N,1);                      % P(rewire|opponentIsCooperator)
    str(:,3) = rand(N,1);                      % P(rewire|opponentIsDefector)
    str(:,4) = s.luceExponentSD*randn(N,1);    % Luce choice exponent
  end
end
