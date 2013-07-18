% Starts as a population of defectors who rewire if the opponent defects,
% but not if they cooperate. Everyone always follows that policy. New mutants
% get a random gameplay strategy. When forming a new link, neighbors are 
% selected with a probability that is a function of their degree. This 
% selection is controlled by a choice exponent. A range of choice exponents 
% is explored. Higher choice exponents (i.e., greater dependence on
% popularity) leads to greater levels of cooperation.
function r = Simulation1()
  s = defaultSettings();            
  r = exploreParameter(s, 'luceMean', [-10,0,10])
  plotTimecourse(r);
end
