function results = exploreParameter(settings,parameter,values)
  for i = 1:length(values)
    if(settings.verbose)
      fprintf('parameter value %d of %d\n', i, length(values))
    end
    % set the parameter value
    if(iscell(values))
      settings = setfield(settings, parameter, values{i});
    else
      settings = setfield(settings, parameter, values(i));
    end   
    % run the simulation       
    parfor j = 1:settings.numRuns
      if(settings.verbose)
        fprintf('  chain %d of %d\n', j, settings.numRuns)
      end
      r = Simulator(settings);
      thisData(:,j) = r.history(:,1);
    end
    results(:,i) = mean(thisData,2);
  end
end