%
%
%
function results = Simulator(s)
  
  % preallocate results matrices
  results.history = zeros(s.numSteps,4);
  
  % initialize the population's interaction graph
  pop.graph = MakeAdjacencyMatrix(s.graphType, s.N, s.graphDensity) > 0;
  
  % initialize the population's strategies
  if(strcmp(s.initial,'random'))
    pop.strategies = generateRandomStrategies(s.N);
  elseif(strcmp(s.initial,'zeroed'))
    pop.strategies(:,1) = zeros(s.N,1);  % probability of cooperating
    pop.strategies(:,2) = zeros(s.N,1);  % P(rewire|opponentDefects)
    pop.strategies(:,3) = zeros(s.N,1);  % P(rewire|opponentCooperates)
    pop.strategies(:,4) = s.luceMean+s.luceSD*randn(s.N,1); % Luce choice exponent
  end
  if(s.isStrategyBinary)
    pop.strategies(:,1) = pop.strategies(:,1) > 0.5;
  end

  for i = 1:s.numSteps
    
    if(rand() < s.pRewireRound) % link update round
      
      % select a random player
      player1 = randi(s.N);
      
      % select one of that player's existing links
      links = find(pop.graph(player1,:));
      if(~isempty(links))
        player2 = links(randi(length(links)));
        
        % the rewiring decision is conditional of player2's behavior
        player2Cooperates = rand() < pop.strategies(player2,1);
        rewire = (rand() < pop.strategies(player1, player2Cooperates+2));
    
        if(rewire)
          % preferential attachment to unlinked players w/ Luce choice exponent
          unlinkedPlayers = find(~pop.graph(player1,:));
          if(~isempty(unlinkedPlayers))
            degree = full(sum(pop.graph, 1));
            L = pop.strategies(player1,4);
            signal = (1+degree(unlinkedPlayers)).^L;
            newLink = unlinkedPlayers(randp(signal./sum(signal)));
    
            % remove the old connection; form the new one
            pop.graph(player1,player2) = 0;
            pop.graph(player1,newLink) = 1;
          end
        end  
      end

    else % strategy update round  
      
      players = shuffle(1:s.N); % student is players(1), teacher is players(2)
      
      if(rand() < s.pMutation)
        pop.strategies(players(1),:) = generateRandomStrategies(1);
      else
        if(strcmp(s.process,'Pairwise'))
          payoffs(1) = playWithNeighbors(players(1)); % student
          payoffs(2) = playWithNeighbors(players(2)); % teacher
        
          % replacement determined by fermi function
          isReplace = rand() < 1/(1+exp(-s.beta*(payoffs(2)-payoffs(1))));
          if(isReplace)
            pop.strategies(players(1),:) = pop.strategies(players(2),:);
          end
        
        elseif(strcmp(s.process,'Moran'))
          for j = 1:s.N
            payoffs(j) = playWithNeighbors(j);
          end
          toReplicate = randp(payoffs./sum(payoffs));
          pop.strategies(players(1),:) = pop.strategies(toReplicate,:);
        end
      end
    end
    % store results
    results.history(i,:) = mean(pop.strategies);
    if(~mod(i,10000) & s.verbose)
      mean(pop.strategies)
    end
  end
  
  % play a game between the player and each of its neighbors
  function payoff = playWithNeighbors(player)
    % compute total cost
    numNeighbors = full(sum(pop.graph(player,:)));
    totalCost = s.cost * numNeighbors * pop.strategies(player,1);
    
    % compute total benefit
    pCooperate = pop.strategies(full(pop.graph(:,player)),1);
    totalBenefit = sum(s.benefit * pCooperate);
    
    % compute final payoff
    payoff = totalBenefit - totalCost;
  end
  
  % TODO: add noise epsilon
  function str = generateRandomStrategies(N)
    str(:,1) = rand(N,1);                      % probability of cooperating
    str(:,2) = rand(N,1);                      % P(rewire|opponentIsDefector)
    str(:,3) = rand(N,1);                      % P(rewire|opponentIsCooperator)
    if(s.isStrategyBinary)
      str(:,1) = str(:,1) > 0.5;
    end
    str(:,4) = s.luceMean+s.luceSD*randn(N,1); % Luce choice exponent
  end
end
