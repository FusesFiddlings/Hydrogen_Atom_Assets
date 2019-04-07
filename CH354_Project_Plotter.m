close all;

plotMargin = 0;
L=[-plotMargin:radialDimension+plotMargin]*stepSize;
W=[-radialDimension-plotMargin:radialDimension+plotMargin]*stepSize;

clearProton=true;
if clearProton
    protonRegion = 2;

    for u=0:protonRegion
        for v=-protonRegion:protonRegion
          if u^2+v^2<protonRegion^2+.1
              parentRadialEnergyMatrix(u+1,v+radialDimension+1)=0;
          end
        end
    end
end

figure('Name','Probability Density');
f=surf(W,L,parentProbabilityMatrix);
set(f,'LineStyle','none');
set(gcf, 'Position',  [100, 100, 1100, 600]);
view(0,90);
drawnow;

%figure('Name','Radial Probability Density');
%f=surf(W,L,parentRadialProbabilityMatrix);
set(f,'LineStyle','none');
set(gcf, 'Position',  [100, 100, 1100, 600]);
view(0,90)
drawnow;

%figure('Name','Energy Density');
%f=surf(W,L,parentEnergyMatrix);
set(f,'LineStyle','none');
set(gcf, 'Position',  [100, 100, 1100, 600]);
view(0,90);
drawnow;

figure('Name','Radial Energy Density');
f=surf(W,L,-1*parentRadialEnergyMatrix);
set(f,'LineStyle','none');
set(gcf, 'Position',  [100, 100, 1100, 600]);
view(0,90);
drawnow;