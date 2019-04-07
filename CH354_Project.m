format compact
format short
clear all;

bohrRadius = 5.29177e-11;
electronCharge = 1.60217662e-19;
coulombsConstant = 8.9875517873681764e9;
vacuumPermittivity = 8.85418782e-12;
reducedMass = 9.1044251e-31;
hBar = 1.0545718e-34;
hydrogenEnergy = reducedMass*electronCharge^4/(8*vacuumPermittivity^2*(hBar*2*pi)^2);
kineticCoefficient = -(hBar^2)/2/reducedMass;
potentialCoefficient = electronCharge^2*coulombsConstant;
 
scale = 5e-10;
radialDimension = 1000;
zDimension = 2*radialDimension;
stepSize = scale/radialDimension;
normalize = true;
 
plotMargin = 0;
L=[-plotMargin:radialDimension+plotMargin]*stepSize;
W=[-radialDimension-plotMargin:radialDimension+plotMargin]*stepSize;
 
zp=2e-10;
charge = 1;
energy = 0;
 
parentEnergy = 1;
energyList = []
 
parameter = 1.01;
eccentricity = 0.09;
 
YBAR = [0:1:5]*stepSize;
YDEV = [3:4:80]*stepSize;
RDEV = [3:4:80]*stepSize;
WE = exp(-4.6);
 
hydrogen = @(x,y,z,p,e)(1./sqrt(pi).*(bohrRadius).^-1.5.*exp(-sqrt(x.^2+y.^2+z.^2)./(p.*bohrRadius.*(1-e.*sin(atan(-y./sqrt(x.^2+z.^2)))))));

hydrogen2 = @(x,y,z)(1./sqrt(pi).*(bohrRadius).^-1.5.*exp(-sqrt(x.^2+y.^2+z.^2)./bohrRadius));
 
gaussian = @(x,y,z,ybar,ydev,rdev)(exp(-(y-ybar).^2/(2*ydev.^2)-(x.^2+z.^2)/(2*rdev.^2)));
 
format long
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
for yBar = YBAR
    for wallExponential = WE
        for yDev = YDEV
            for rDev = RDEV
            wFunc=zeros(radialDimension+3,zDimension+3,3);
            
            radialProbabilityMatrix = zeros(radialDimension+3,zDimension+3);
            probabilityMatrix = radialProbabilityMatrix;
 
            for u = [-1:radialDimension+1]
                for v = [-radialDimension-1:radialDimension+1]
                    for w = [-1:1] 
                        %if ~(u==0 && v==0 && w==0)
                         %   hVal=hydrogen(u*stepSize,v*stepSize,w*stepSize,parameter,eccentricity);
                        %else
                            %hVal=hydrogen2(u*stepSize,v*stepSize,w*stepSize);
                        %end
                        
                        if v*stepSize+zp <= 0 && charge > 0
                            wFunc(u+2,v+radialDimension+2,w+2) = 0;
                        else 
                            %wFunc(u+2,v+radialDimension+2,w+2) = (hVal*stepSize^1.5)*(1-exp(-wallExponential*(v*stepSize+zp)/stepSize));
                            wFunc(u+2,v+radialDimension+2,w+2) = (hVal*stepSize^1.5)*gaussian(u*stepSize,v*stepSize,w*stepSize,yBar,yDev,rDev)*(1-exp(-wallExponential*(v*stepSize+zp)/stepSize));
                        end
                        if(w==0 && u>-1)
                            radialProbabilityMatrix(u+2,v+radialDimension+2)=wFunc(u+2,v+radialDimension+2,w+2)^2*2*pi*u;
                            probabilityMatrix(u+2,v+radialDimension+2)=wFunc(u+2,v+radialDimension+2,w+2)^2;
                        end 
                    end
                end
            end
            probability=sum(sum(radialProbabilityMatrix));
 
            if normalize
                wFunc=wFunc/sqrt(probability);
                for u = [-1:radialDimension+1]
                    for v = [-radialDimension-1:radialDimension+1]
                        radialProbabilityMatrix(u+2,v+radialDimension+2)=wFunc(u+2,v+radialDimension+2,2)^2*2*pi*u;
                    end
                end
                renormalizedProbability=sum(sum(radialProbabilityMatrix));
            end
            
            
            radialEnergyMatrix = zeros(radialDimension+1,zDimension+1);
            energyMatrix = radialEnergyMatrix;
            
            for u = [0:radialDimension]
                for v = [-radialDimension:radialDimension]
                    currentwFunc = wFunc(u+2,v+radialDimension+2,2);
                    wdXderiv = (wFunc(u+1,v+radialDimension+2,2)-2*currentwFunc+wFunc(u+3,v+radialDimension+2,2))/stepSize^2;
                    wdYderiv = (wFunc(u+2,v+radialDimension+1,2)-2*currentwFunc+wFunc(u+2,v+radialDimension+3,2))/stepSize^2;
                    wdZderiv = (wFunc(u+2,v+radialDimension+2,1)-2*currentwFunc+wFunc(u+2,v+radialDimension+2,3))/stepSize^2;
 
                    if (u>0 || u<0 || v>0 ||v<0) && zp+v*stepSize>0
                        kineticEnergy = kineticCoefficient*(wdXderiv+wdYderiv+wdZderiv);
                        potentialEnergy = potentialCoefficient*(-1/(stepSize*sqrt(u^2+v^2))+charge/(zp+v*stepSize));
                        radialEnergyMatrix(u+1,v+radialDimension+1) = 2*pi*u*currentwFunc*(kineticEnergy+currentwFunc*potentialEnergy);
                        energyMatrix(u+1,v+radialDimension+1) = currentwFunc*(kineticEnergy+currentwFunc*potentialEnergy);
                    end
                end
            end
 
            energy = sum(sum(radialEnergyMatrix));
            %energyList = [energyList energy];
            
            if energy < parentEnergy
                parentRadialProbabilityMatrix = radialProbabilityMatrix(2:2+radialDimension,2:2+zDimension);
                parentProbabilityMatrix = probabilityMatrix(2:2+radialDimension,2:2+zDimension);
                parentRadialEnergyMatrix=radialEnergyMatrix;
                parentEnergyMatrix=energyMatrix;
                parentEnergy = energy
                bestCombo = [log(wallExponential) yBar/stepSize yDev/stepSize rDev/stepSize]
 
                run('CH354_Project_Plotter.m')
            end
        end
    end
    end
end
beep on
beep


