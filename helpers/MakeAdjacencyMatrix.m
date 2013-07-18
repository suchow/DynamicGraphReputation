% Returns the transition matrix for the given graph type
function [w,varargout] = MakeAdjacencyMatrix(graphType,N,varargin)
  
  switch graphType

  % Complete graph with identically-weighted edges, as in the Moran Process
  case 'Complete'
    if(N == 1)
      w = 1;
    else
      w = ~diag(1:N)./(N-1);
    end  
  
  % A complete bipartite graph
  case 'Complete bipartite'
    n = varargin{1};
    m = varargin{2};
    if(n+m ~= N)
      error('In a bipartite graph, the two sides must sum to N.')
    end
    w = sparse(N,N);
    w(1:n,n+1:end) = 1;
    w((n+1):N,1:n) = 1;
    w = rownormalize(w);
  
  
  % The G(n, p) Erdős–Rényi model proposed in Edgar Gilbert's 1959 "Random
  % graphs", except with no isolated vertices. Defaults to p = 0.5 unless
  % specified as an argument to MakeAdjacencyMatrix.
  case 'Erdos-Renyi'
    if(length(varargin) < 1)
      p = 0.5;
    else
      p = varargin{1};
    end
    w = sparse(N,N);
    while(~all(sum(w,2) ~= 0)) % resample all empty rows
      emptyRows = (sum(w,2) == 0);
      w(emptyRows,:) = (rand(sum(emptyRows),N) < p);
      w2 = w'; % make symmetric step 1
      t = triu(true(size(w2)),1);
      w(t) = w2(t);
      w = w - diag(diag(w)); % zero out diagonal    
    end
    w = rownormalize(w);

  % Square lattice connected to its von Neumann neighborhood
  case 'Square lattice 4'
    r = floor(sqrt(N));
    c = floor(sqrt(N));
    diag1 = repmat([ones(c-1,1); 0],r,1);  
    diag1 = diag1(1:end-1);
    diag2 = ones(c*(r-1),1);
    w = diag(diag1,1) + diag(diag2,c);
    w = w+w.';  
    w = sparse(rownormalize(w));  


  % Square lattice connected to its Moore neighborhood
  case 'Square lattice 8'
    r = floor(sqrt(N));
    c = floor(sqrt(N));
    diag1 = repmat([ones(c-1,1); 0],r,1);
    diag1 = diag1(1:end-1);            
    diag2 = [0; diag1(1:(c*(r-1)))]; 
    diag3 = ones(c*(r-1),1);         
    diag4 = diag2(2:end-1);       
    w = diag(diag1,1) + diag(diag2,c-1) + ...
        diag(diag3,c) + diag(diag4,c+1);
    w = w+w.';                    
    w = sparse(rownormalize(w));  
          

  % A path from vertex 1 to N
  case 'Path'
    w = sparse(N,N);
    w = addlinks(w, 1:N, min(2:(N+1),N));

  % A directed cycle created by having the snake eat its tail
  case 'Cycle'
    w = sparse(N,N);
    w = linktoneighbors(w,1);


  case 'Star'
    w = sparse(N,N);
    w(1,2:end) = true;
    w(2:end,1) = true;
    w = rownormalize(w);


  case 'Burst'
    w = sparse(eye(N));
    w(1,1) = false;
    w(1,2:end) = true;
    w = rownormalize(w);

  % All nodes connect to a central hub
  case 'Integrator'
    w = sparse(N,N);
    w(:,1) = true;


  % Uses the Barabasi-Albert model to construct a scale free network with
  % scaling exponent 3. We seed it with m_0 = 2 nodes, and an edge from 
  % node 1 to node 2 and back. Each step connects a new node to two existing
  % nodes with probability proprtional to their degree
  case 'Scale free'
    w = sparse(N,N);
    w(1,2) = true;
    w(2,1) = true;
    w(3,1) = true;
    w(1,3) = true;
    w(2,3) = true;
    w(3,2) = true;
    for i = 4:N
      p_i = rownormalize(sum(full(w)')); % count edges
      n1 = randp(p_i);          % pick first node
      p_i(n1) = 0;                       % disallow repeat (w/o replacement)
      p_i = rownormalize(p_i);           % renormalize
      n2 = randp(p_i);          % pick second node
      w(i,n1) = true;                    % link up the chosen nodes
      w(n1,i) = true;
      w(i,n2) = true;
      w(n2,i) = true;
    end
    w = rownormalize(w);
  
  
  % Barabasi-Albert Model A (growth, but no preferential attachment)
  case 'Model A'
    w = sparse(N,N);
    w(1,2) = true; % initialize with 3 nodes
    w(2,1) = true;
    w(3,1) = true;
    w(1,3) = true;
    for i = 4:N
      x = shuffle(1:(i-1));
      w(i,x(1)) = true;  % link up the chosen nodes
      w(x(1),i) = true;
      w(i,x(2)) = true;
      w(x(2),i) = true;
    end
    w = rownormalize(w);
  
  
  case 'Friendship'
    if(~mod(N,2)) % if odd
      error('Friendship graph must have an odd number of vertices');
    end
    numMills = (N - 1) / 2;
    w = sparse(N,N);
    w(1,2:end) = true;
    for i = 1:numMills
      w(2*i,2*i+1) = true;
      w(2*i+1,2*i) = true;
      w(2*i,1) = true;
      w(1,2*i) = true;
      w(2*i+1,1) = true;
      w(1,2*i+1) = true;
    end
    w = rownormalize(w);
    
    
  % From Mertzios et al. (2012)  
  case 'Clique wheel'
    w = sparse(N,N);
    n = N/2;
    w(1:n,1:n) = linktoneighbors(w(1:n,1:n),[1,-1]);
    w((n+1):N,(n+1):N) = linktoneighbors(w((n+1):N,(n+1):N),[1,-1]);
    for i = 1:n % link inner to outer wheel
      w(i,n+i) = true;
      w(n+i,i) = true;
    end
    w = rownormalize(w);
    
    
  case 'Konigsberg'
    w = sparse(4,4);
    w(1,2) = 2;
    w(2,1) = 2;
    w(1,3) = 2;
    w(3,1) = 2;
    w(1,4) = 1;
    w(4,1) = 1;
    w(3,4) = 1;
    w(4,3) = 1;
    w(2,4) = 1;
    w(4,2) = 1;
    w = rownormalize(w);


  % Watts and Strogatz model
  % note: this algorithm uses a different order for rewiring.
  % ** extend to arbitrary k
  case 'Watts-Strogatz'
    if(nargin < 3)
      p = 0.5;
    else
      p = varargin{1}/2;
    end
    w = sparse(N,N);
    w = linktoneighbors(w,[1,2,-1,-2]);
    % random rewiring
    [i,j] = find(w);
    for k = 1:length(i)
      if((rand < p) && (j(k) > i(k)))
        % break the old link between i and j
        w(i(k),j(k)) = false;
        w(j(k),i(k)) = false;
        % choose a new target
        newTarget = randi(N);
        % create link from i to target
        w(i(k),newTarget) = true;
        w(newTarget,i(k)) = true;
      end
    end
    w = rownormalize(w);
      
  
  % scale free model with gaussian fitness, as in Bianconi & Barabasi
  case 'Scale free (GF)'
    w = sparse(N,N);
    mu = 1; % mean log fitness
    sd = 2; % sd of log fitness
    f = lognrnd(mu,sd,1,N);
    w(1,2) = true;
    w(2,1) = true;
    w(3,1) = true;
    w(1,3) = true;
    for i = 4:N
      p_i = rownormalize(f .* rownormalize(sum(full(w)'))); % count edges
      n1 = randp(p_i);          % pick first node
      p_i(n1) = 0;                       % disallow repeat (w/o replacement)
      p_i = rownormalize(p_i);           % renormalize
      n2 = randp(p_i);          % pick second node
  
      w(i,n1) = true;                    % link up the chosen nodes
      w(n1,i) = true;
      w(i,n2) = true;
      w(n2,i) = true;
    end
    w = rownormalize(w);
    
    
  % Randomly diluted lattice, from Dhar et al. (1987)
  % removes nodes indepdendently with probability p
  case 'Diluted square lattice 4'
    if(nargin < 3)
      p = 0.75;
    else
      p = varargin{1};
    end
    w =  MakeAdjacencyMatrix('Square lattice 4', N) > 0;
    w = dilute(w,p);
    w = rownormalize(w);
  
  
  case 'Diluted square lattice 8'
    if(nargin < 3)
      p = 0.75;
    else
      p = varargin{1};
    end
    w =  MakeAdjacencyMatrix('Square lattice 8', N) > 0;
    w = dilute(w,p);
    w = rownormalize(w);
    
    
  % ADD: random grid from Nowak's 1994 paper 'More spatial games'
  % case 'Random grid'
    
  % The Bethe lattice, also called the Cayley tree, from Bethe (1935)
  case 'Bethe lattice'
    if(nargin < 3)
      k = 3;
    else
      k = varargin{1};
    end
    w = sparse(N,N);
    i = 1;
    newNode = 3;
    w(1,2) = true; % link up the first node
    w(2,1) = true;
    while(newNode <= N)
      for j = 1:(k-1)
        w(i,newNode) = true;
        w(newNode,1) = true;
        newNode = newNode + 1;
      end
      i = i + 1;
    end
    w = rownormalize(w);
    
  
  % Random regular graph (each node has the same number of edges)
  % Uses the method of Kim & Vu (2006)
  case 'Random regular'
    if(nargin < 3) % the degree
      d = 3; 
    else
      d = varargin{1}; 
    end
    w = sparse(N,N);
    tmp = repmat(1:N,d,1);
    U = tmp(:)';
    Uw = sparse(N,N); % keep track of who we have linked
    pairings = [];
    while(~isempty(U))
      ord = shuffle(1:length(U));
      i = U(ord(1));
      j = U(ord(2));
      % check if the pair is suitable:
      isLoop = (i==j);
      isParallel = Uw(i,j);
      suitable = ~isLoop & ~isParallel;
      if(suitable)
        Uw(i,j) = true;
        Uw(j,i) = true;
        U([ord(1) ord(2)]) = [];
      end
    end
     % if Uw is regular, output it
    isRegular = all(sum(Uw) == d) & all(sum(Uw') == d);
    if(isRegular)
      w = Uw;
    else
      w = MakeAdjacencyMatrix('Random regular', N, d);
    end    
    w = rownormalize(w);
    

  % Every node is an island
  case 'Isolated'
    w = speye(N,N);
  
  % ADD: pseudofractal scale-free web
  % http://arxiv.org/pdf/cond-mat/0112143v1.pdf
  
  % Constructs a graph using the Luce choice axiom
  % L can be an XRP or Inf.
  case 'Luce'    
    w = sparse(N,N);
    m = varargin{1};
    L = varargin{2};
    if(nargin < 5)
      v_function = @(w)(rownormalize(full(sum(w))));
    else
      v_function = varargin{3};
    end
    %v_function = @(w)(pagerank(w')');
    if(isnumeric(L))
      L = @()(L);
    end
    % start with a complete graph with m nodes
    w(1:(m+1),1:(m+1)) = MakeAdjacencyMatrix('Complete', m+1);
    personalL = nan(1,N);  
    for i = (m+2):N
      if(N>10000 && ~mod(i,100))
        fprintf('%d\n',i)
      end
      used = i;
      personalL(i) = L(); % each individual gets assigned a value of L
      for j = 1:m
        v = v_function(w(1:(i-1),1:(i-1))); % find value of each node
        p = zeros(1,N);
        if(isinf(personalL(i)))
          v_tmp = v;
          if(personalL(i) > 0) % +Inf case
            v_tmp(used) = 0;
            ms = (v_tmp==max(v_tmp));
          else % -Inf case
            v_tmp(~v_tmp) = Inf;
            v_tmp(used) = Inf;
            ms = (v_tmp==min(v_tmp));
          end
          p(ms) = 1./sum(ms);
        else % -Inf < L < Inf case
          p(v>0) = v(v>0).^personalL(i); % raise each to the Lth power
          p(used) = 0;
          p = p./sum(p);
        end
        x = randp(p); % pick a node
        used = [used x];
        w(i,x) = true;
        w(x,i) = true;
      end
    end
    w = rownormalize(w);
    varargout{1} = personalL;
    
  
  case 'Visual'
    r = floor(sqrt(N));
    c = floor(sqrt(N));
    N = r*c;
    %h = varargin{1};
    h = 3;
    w = zeros(N,N);
    for i = 1:r
      for j = 1:c
        [x,y] = meshgrid(1:r,1:c);
        d = ((x-i).^2 + (y-j).^2).^0.5; % compute distance
        d2 = exp(1).^-(d.^2/(2*h.^2)); % compute strength of connection
        w(sub2ind([r,c],i,j),:) = d2(:)';
      end
    end
    w(diag(1:N)>0) = 0;
    w = rownormalize(w);    
  
  otherwise
    error('This graph type is not supported.')
    
  end
end

function w = addlinks(w,from,to)
  w(sub2ind(size(w),from,to)) = true;
end

function w = linktoneighbors(w,t)
  N = length(w);
  for i = 1:length(t)
    w = addlinks(w,1:N,1+mod(t(i)+(0:(N-1)),N));
  end
end

% remove links with probably p, and relink isolated nodes to self
function w = dilute(w,p)
  kill = rand(1,length(w)) < 1-p;
  w(kill,:) = [];
  w(:,kill) = [];
  % link isolated nodes to self
  islands = ~sum(w);
  w = diag(islands) | w;
end
