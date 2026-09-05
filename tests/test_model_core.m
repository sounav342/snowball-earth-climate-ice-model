function tests = test_model_core
tests = functiontests(localfunctions);
end

function testFractionalAlbedoLimits(testCase)
par = struct('nj',3,'alpha_i',0.6,'alpha_o',0.1);
R = [0;0.5;1];
alpha = surface_albedo_fractional(R,par);
verifyEqual(testCase,alpha,[0.1;0.35;0.6],'AbsTol',1e-14);
end

function testWeightedAlbedoMatchesWeightedFlux(testCase)
par = struct('nj',5,'alpha_i',0.6,'alpha_o',0.1);
R = linspace(0,1,5).';
SWdown = [100;150;200;250;300];
alpha = surface_albedo_fractional(R,par);
fromAlbedo = SWdown .* (1-alpha);
fromFluxes = R .* SWdown .* (1-par.alpha_i) ...
    + (1-R) .* SWdown .* (1-par.alpha_o);
verifyEqual(testCase,fromAlbedo,fromFluxes,'AbsTol',1e-12);
end

function testOutOfRangeIceFractionFails(testCase)
par = struct('nj',2);
verifyError(testCase,@() validate_fractional_ice_cover([0;1.01],par,'test'), ...
    'validate_fractional_ice_cover:outOfRange');
end
