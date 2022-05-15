%% Calibration of the Ocean spectrometer

%This section shows the calibration  of the ocean spectrometer using the
%calibration source. 
Expected1 = [365.015
404.656
407.783
435.833
546.074
576.96
579.066
696.543
706.722
727.294
738.393
750.387
763.511
772.376
794.818
800.616
811.531
826.452
852.144
866.794
912.297
922.45];

Recorded1 = [367.11
406.46
409.60
438.41
549.29
578.51
580.62
697.97
708.17
728.82
739.84
751.76
764.92
773.91
796.22
802.32
812.81
827.80
853.51
868.20
913.49
923.72];

Difference1 = Recorded1 - Expected1;
%This equation gives the delta between the expected values and the
%recorded values so that we can calculate the standard deviation.

Average = mean(Difference1);
Standard_derivation = std(Difference1);

Standard_error = Standard_derivation/sqrt(length(Difference1));

sserrorplot = 0.1*ones(length(Difference1),1);

hold on
errorbar(Expected1,Difference1,sserrorplot,'.b','MarkerSize',8,'CapSize',4)
xlabel('Wavelength (nm)','FontSize', 20,'Interpreter','latex')
ylabel('Difference in wavelength (nm)')
ylim([0 4])
xlim([350 950])
set(gca,'FontSize',16)
legend('off')
hold off

% This is the associated error on the average
%% Calibration of the Bentham monochromator

clc

% Inputting data manually recorded
Expected2 = [
400
425
450
475
500
525
550
575
600
625
650
675
700
725
750
775
800
825
850
875
900
925
950
975
1000];

Recorded2 = [416.51
441.12
465.80
490.33
514.93
539.97
564.64
589.53
614.23
639.14
663.84
688.71
713.54
738.33
763.24
788.07
812.81
837.62
862.48
887.02
911.57
936.12
960.96
985.74
1010.74];

Adjusted2 = Recorded2 - Average;

Difference2 = Adjusted2 - Expected2;

standard_deviation2 = std(Difference2);
weights = 1 / standard_deviation2.^2; 
weightsplot = weights*ones(length(Difference2),1);

%Using the standard error approach to combine the standard error on the

serror = sqrt(0.1^2 + Standard_error^2);
serrorplot = serror * ones(length(Difference2),1);
clf
hold on
errorbar(Expected2,Difference2,serrorplot,'*b','MarkerSize',10)
h = plot(fittedmodel1);
xlabel('Wavelength (nm)','FontSize', 20)
ylabel('Difference in wavelength (nm)')
xlim([350 1050])
ylim([8 16])
set(gca,'FontSize',16)
set(h,'LineWidth',2)
legend('off')
hold off

%Pulling coeffecients and errors
Bound = confint(fittedmodel1);
Gradlow = Bound(1,2);
Gradup = Bound(2,2);
Intlow = Bound(1,1);
Intup = Bound(2,1);

Grad_err = (Gradup - Gradlow)/2*1.96
Int_err = (Intup - Intlow)/2*1.96





