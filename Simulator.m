function results = Simulator(s)
  
  % derivative parameters
  isRewireRound = rand(1,s.numRounds) < s.pRewireRound;

  % create the initial pop
  pop.graph = MakeAdjacencyMatrix(s.graphType, s.N);
  pop.strategies = generateRandomStrategies(s.N); %
  
  for i = 1:s.numRounds
    
    % rewire round
    if(rand < s.pRewireRound)
      % pick a random connection, player 1 makes the choice
      players = shuffle(1:s.N);
      player1strategy = pop.strategies(players(1),:);
      player2strategy = pop.strategies(players(2),:);
      
      rewire = rand() < player1strategy(player2strategy(1)+1);
      if(rewire)
        % preferential attachment to unlinked players w/ Luce choice exponent
        unlinkedPlayers = find(~pop.graph(players(1),:));    
        degree = full(sum(pop.graph > 0, 2));
        L = player1strategy(4);
        p = (degree(unlinkedPlayers).^L)/sum(degree(unlinkedPlayers).^L);
        newConnection = unlinkedPlayers(randp(p));
        
        % remove the old connection; form the new one
        pop.graph(players(1),players(2)) = 0;
        pop.graph(players(1),newConnection) = 1;
        pop.graph = rownormalize(pop.graph);
      end

    % playing round  
    else
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
        if(strcmp(s.process,'Pairwise'))
          % simulate a pairwise comparison process: choose a pair of players,
          % assign one the role of student, the other the role of teacher.
          % the student takes on the teacher's strategy with probability
          % proportional to payoff. implicit in this model is a replacement
          % network that is complete.
          players = shuffle(1:s.N);
          payoffs(1) = playWithNeighbors(players(1)); % student
          payoffs(2) = playWithNeighbors(players(2)); % teacher
        
          isReplace = rand() < payoffs(2)/sum(payoffs(1:2));
          if(isReplace)
            pop.strategies(players(1),:) = pop.strategies(players(2),:);
          end
        
        elseif(strcmp(s.process,'Moran'))
          for j = 1:s.N
            payoffs(j) = playWithNeighbors(j);
          end
          toReplicate = randp(payoffs./sum(payoffs));
          pop.strategies(randi(s.N),:) = pop.strategies(toReplicate,:);
        end
      end
    end
    % store results
    results.history(i) = pop;
  end
  
  % play a game between player and each of its neighbors
  function payoff = playWithNeighbors(player)
    neighbors = find(pop.graph(player,:));
    for k = 1:length(neighbors)
      payoff(k) = s.payoffMatrix(pop.strategies(player,1),...
                                 pop.strategies(neighbors(k),1));
    end
    payoff = sum(payoff);
  end
  
  % TODO: add noise epsilon
  function str = generateRandomStrategies(N)
    str(:,1) = 1 + (rand(N,1) > s.pCooperator);% allC = 1, allD = 2
    str(:,2) = rand(N,1);                      % P(rewire|opponentIsCooperator)
    str(:,3) = rand(N,1);                      % P(rewire|opponentIsDefector)
    str(:,4) = s.luceMean+s.luceSD*randn(N,1);    % Luce choice exponent
  end
end
