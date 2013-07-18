% Consider the effect of different cost/benefit ratios
function r = Simulation4()
  s = defaultSettings();            
  s.luceMean = 10;
  s.cost = 1;
  r = exploreParameter(s, 'benefit', [2,10,20,100,200,1000])
  plotTimecourse(r)
end
