function results = Simulation1

	% settings
	s.N = 100;                             % pop size
	s.graphType = 'Erdos-Renyi';           % initial network structure
	s.numRounds = 2000;                    % total number of rounds
	s.payoffMatrix = [3, 0; 5, 1];         % payoff matrix for the game
	s.pCooperator = 0.5;                   % initial probability of playing allC
	s.pRewireRound = 0.1;                  % probability of a rewire round
	s.pMutation = 0.1;                     % probability of mutating on play round
	s.luceExponentSD = 2;                  % standard deviation of Luce choice exp
	
	results = Simulator(s)

end