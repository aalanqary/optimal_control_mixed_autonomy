function fc = simplified_fuel_model(v,a,vehicle)

% This function implements the formulae for the simplified
% fuel consumption models of various vehicles.
% (C) 2020/12/23 by Rabie Ramadan, Benjamin Seibold,
% Kenneth Butts, Joy Carpio , and Sulaiman Almatrudi

% output1 = fc = fuel consumption (in grams/sec)
% output2 = P = power (in KW)
% output3 = flag = {0 if no warning, 1 if input is outside
%                   model region of feasibility, 2 if velocity
%                   is negative}
%
% NOTE: when velocity is negative, it will be treated as a zero
% for fuel consumption calculation purposes. When input is outside
% model region of feasibility, fc will be assigned a value based
% on the fuel consumption polynomial, which extends gracefully
% into the region of infeasiblity.
%
% input1 = v = vehicle velocity (in m/s)
% input2 = a = vehicle acceleration (in m/s^2)
% input3 = vehicle = vehicle name (one of the following strings)
% {'midBase','RAV4','midSUV','Pickup','Class3PND','Class8Tractor'}
% Vehicle descriptions:
% midBase       <-- This is a mid size compact sedan with mass 1743 kg
% RAV4          <-- This is the 2019 Toyota Rav4 with mass 1717 kg
% midSUV        <-- This is a mid size sports SUV with mass 1897 kg
% Pickup        <-- This is a light-duty pickup truck with mass 2173 kg
% Class3PND     <-- This is a class 3 truck used for pickup and delievery
%                   with mass 5943 kg
% Class8Tractor <-- This is a class 8 tractor truck with manual transmission
%                   It also has a trailer and total mass 25104 kg
%
% input1 and input2 can be matrices or vectors of the same size
% the outputs fc, P and flag will have the same size matrix/vector
% as the inputs.

load([vehicle '_coeffs.mat'])       % load vehicle coefficients
C2 = double(C2);
q0 = double(q0);
gs2KW = 37.42;                      % grams/sec to KW conversion factor
v_max_fit = 40;                     % maximum speed considered in fitting

ma = b1 * (v/v_max_fit).^b2 .* ...  % calculate max feasible
    (1-v/v_max_fit).^b3 +b4*v + b5; % a for each input v

v = max(v,0);                       % treat negative velocities as zero
fc = C0 + C1*v + C2*v.^2 + ...      %
    C3*v.^3 + p0*a + ...            % 
    p1*a.*v + p2*a.*v.^2;% + ...      % polynomial fuel consumption formula 
%     q0*max(a,0).^2 + ...            %
%     q1*max(a,0).^2.*v;              %
fc = max(fc, beta0);                % assign min fc when polynomial is below the min




