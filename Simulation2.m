% Starts as a population of defectors, who rewire depending on the given
% policy: (1) break links with both defectors and cooperators, (2) break links 
% with cooperators, but never defectors, (3) break links with defectors, but 
% never cooperators, or (4) never break links.
function r = Simulation2()
  s = defaultSettings();            
  s.luceMean = 10;
  r = exploreParameter(s, 'fixedRewiring', {[1 0], [0 1], [1 1], [0 0]})
  plotTimecourse(r)
end
