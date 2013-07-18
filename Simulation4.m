% Consider the effect of different rates of rewiring.
function r = Simulation4()
  s = defaultSettings();            
  s.luceMean = 1;
  s.numSteps = 10000;
  r = exploreParameter(s, 'pRewireRound', [0.10, 0.25, 0.50, 0.75, 0.90])
  plotTimecourse(r)
end
