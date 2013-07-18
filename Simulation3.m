% Just like Simulation #2, but instead of trying only 4 rewiring strategies 
% and plotting the full timecourse of the emergence of cooperation, try all
% possible conditional rewiring strategies and plot only the steady-state
% level of cooperation. 
function r = Simulation3()
  s = defaultSettings();            
  s.luceMean = 10;
  s.numSteps = 1e5;
  s.numRuns = 2;
  s.rewireSteps = [0:0.1:1]
  s.lastStepsToAverage = 1e4;
  for i = 1:length(s.rewireSteps)
    for j = 1:length(s.rewireSteps)
      rwD = s.rewireSteps(i);
      rwC = s.rewireSteps(j);
      rTmp = exploreParameter(s, 'fixedRewiring', {[rwD rwC]});
      r(i,j) = median(rTmp((s.numSteps-s.lastStepsToAverage):end));
      figure(10)
      imagesc(r)
      colormap(gray)
      colorbar
      drawnow
  end
  r
end
