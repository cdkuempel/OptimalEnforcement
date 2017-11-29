function Create_diluting_plot()

% Specify which levels of protection to search over
%In this example we search from 0 to 85% protection in 50 increments 
Protection = [linspace(0,0.1,25) linspace(0.1,0.85,25)]; 

% This is the resolution of the search algorithm. This should be >= 1000.
RES = 1000; 

% Go through levels of protection, one-by-one
for p = 1:length(Protection)
   % For each level of protection, identify the harvested population equilibrium
   TB(p) = submodel(Protection(p),RES);
end

% Plot the results
figure(1), clf, hold on, box on; FS = 15;
set(gca, 'fontsize', 14,'ytick',[100:10:200])
plot(Protection,100*TB./min(TB),'k','linewidth',3)
xlabel('Proportion of patches protected','fontsize',FS+4)
ylabel('Relative metapopulation abundance','fontsize',FS+4)
xlim([-0.04 0.85])
ylim(100*[0.95 1.55])
