function tests = CRLBTest
% overall CRLB Test Function
tests = functiontests(localfunctions);
end


function setupOnce(testCase)
% add CRLB functions to path
% first add path
testCase.TestData.origPath = pwd;
testCase.TestData.tmpFolder = '../CommonFunctions';  % common functions are in folder above
cd(testCase.TestData.tmpFolder)
end

function teardownOnce(testCase)
cd(testCase.TestData.origPath)
end

function setup(testCase)
% initialise code we will test
testCase.TestData.CRs = CramerRaoFunctions;
testCase.TestData.cv = 1;
testCase.TestData.wavelength = 309.9604;
testCase.TestData.T = 293.15;
testCase.TestData.n0 = ...
    testCase.TestData.CRs.n_m(testCase.TestData.wavelength,testCase.TestData.T,testCase.TestData.cv);
end

function teardown(testcase)
clearvars testCase.TestData.CRs
clearvars testCase.TestData.cv
clearvars testCase.TestData.wavelength
clearvars testCase.TestData.T
end

function testJthetathetaf(testCase)
[~, ~, A, B, C, ~] = testCase.TestData.CRs.InstrResp(1,1.3,testCase.TestData.n0);
theta = pi./4;
phi = pi./2;
time = 1;
Beta = ([3; 3; 3; 3]);
I = ([25; 25; 25; 25]);
verifyEqual(testCase, ...
    testCase.TestData.CRs.Jthetathetaf(theta, phi, time, Beta, I, A, B, C, 4),...
    2.2548,'AbsTol', 1e-4, 'Jthetatheta not behaving as expected')
end

function testJthetathetaNB3f(testCase)
[~, ~, A, B, C, ~] = testCase.TestData.CRs.InstrResp(1,1.3,testCase.TestData.n0);
theta = pi./4;
phi = pi./2;
time = 1;
I = ([25; 25; 25; 25]);
verifyEqual(testCase, ...
    testCase.TestData.CRs.JthetathetafNB(theta, phi, time, I, A, B, C, 3),...
    16.3171,'AbsTol', 1e-4, 'JthetathetafNB not behaving as expected, 3 Detector Case')
end

function testJthetathetaNB4f(testCase)
[~, ~, A, B, C, ~] = testCase.TestData.CRs.InstrResp(1,1.3,testCase.TestData.n0);
theta = pi./4;
phi = pi./2;
time = 1;
I = ([25; 25; 25; 25]);
verifyEqual(testCase, ...
    testCase.TestData.CRs.JthetathetafNB(theta, phi, time, I, A, B, C, 4),...
    16.3171,'AbsTol', 1e-4, 'JthetathetafNB not behaving as expected, 4 Detector Case')
end

function testJthetaphif(testCase)
[~, ~, A, B, C, ~] = testCase.TestData.CRs.InstrResp(1,1.3,testCase.TestData.n0);
theta = 1.25;
phi = 0.5;
time = 1;
I = ([25; 25; 25; 25]);
Beta = ([3; 3; 3; 3]);
verifyEqual(testCase, ...
    testCase.TestData.CRs.Jthetaphif(theta, phi, time, Beta, I, A, B, C, 4),...
    0.0339,'AbsTol', 1e-4, 'Jthetaphif not behaving as expected, 4 Detector Case')
end

function testJthetaphiNBf(testCase)
[~, ~, A, B, C, ~] = testCase.TestData.CRs.InstrResp(1,1.3,testCase.TestData.n0);
theta = pi./4;
phi = pi./3;
time = 1;
I = ([25; 25; 25; 25]);
verifyEqual(testCase, ...
    testCase.TestData.CRs.JthetaphifNB(theta, phi, time, I, A, B, C, 3),...
    4.80397,'AbsTol', 1e-5, 'JthetaphifNB not behaving as expected, 3 Detector Case')
verifyEqual(testCase, ...
    testCase.TestData.CRs.JthetaphifNB(theta, phi, time, I, A, B, C, 4),...
    -1.97489,'AbsTol', 1e-5, 'JthetaphifNB not behaving as expected, 4 Detector Case')
end


function testJphiphif(testCase)
[~, ~, A, B, C, ~] = testCase.TestData.CRs.InstrResp(1,1.3,testCase.TestData.n0);
theta = pi./4;
phi = pi./2;
time = 1;
I = ([25; 25; 25; 25]);
Beta = ([3; 3; 3; 3]);
verifyEqual(testCase, ...
    testCase.TestData.CRs.Jphiphif(theta, phi, time, Beta, I, A, B, C, 4),...
    4.51353,'AbsTol', 1e-5, 'Jphiphif not behaving as expected, 4 Detector Case')
end

function testJphiphiNB4f(testCase)
[~, ~, A, B, C, ~] = testCase.TestData.CRs.InstrResp(1,1.3,testCase.TestData.n0);
theta = pi./4;
phi = pi./2;
time = 1;
I = ([25; 25; 25; 25]);
verifyEqual(testCase, ...
    testCase.TestData.CRs.JphiphifNB(theta, phi, time, I, A, B, C, 4),...
    20.3109,'AbsTol', 1e-4, 'JphiphifNB not behaving as expected, 4 Detector Case')
end

function testJphiphiNB3f(testCase)
[~, ~, A, B, C, ~] = testCase.TestData.CRs.InstrResp(1,1.3,testCase.TestData.n0);
theta = pi./4;
phi = pi./2;
time = 1;
I = ([25; 25; 25; 25]);
verifyEqual(testCase, ...
    testCase.TestData.CRs.JphiphifNB(theta, phi, time, I, A, B, C, 3),...
    10.1554,'AbsTol', 1e-4, 'JphiphifNB not behaving as expected, 3 Detector Case')
end

function testQthetatheta1s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./2;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
Beta = 1e4;
I = 500;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qthetatheta1(theta, phi, Beta, I, A, B, C, H, a11, a13, a33),...
    0.1166,'AbsTol', 1e-4, 'Qthetatheta1 not behaving as expected')
end

function testQthetatheta1NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./2;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
I = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qthetatheta1NB(theta, phi, I, A, B, C, H, a11, a13, a33),...
    2.33325,'AbsTol', 1e-5, 'Qthetatheta1NB not behaving as expected')
end

function testMthetatheta1s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./2;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
Beta = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mthetatheta1(theta, phi, Beta, A, B, C, H, NF, a11, a13, a33),...
    0.00893141,'AbsTol', 1e-8, 'Mthetatheta1 not behaving as expected')
end

function testMthetatheta1NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./2;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mthetatheta1NB(theta, phi, A, B, C, H, NF, a11, a13, a33),...
    0.00892771,'AbsTol', 1e-8, 'Mthetatheta1NB not behaving as expected')
end

function testQthetatheta2s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
Beta = 1e4;
I = 5000;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qthetatheta2(theta, phi, Beta, I, A, B, C, H, a11, a13, a33),...
    0.218698,'AbsTol', 1e-6, 'Qthetatheta2 not behaving as expected')
end

function testQthetatheta2NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
I = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qthetatheta2NB(theta, phi, I, A, B, C, H, a11, a13, a33),...
    0.437484,'AbsTol', 1e-6, 'Qthetatheta2NB not behaving as expected')
end

function testMthetatheta2s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./2;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
Beta = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mthetatheta2(theta, phi, Beta, A, B, C, H, NF, a11, a13, a33),...
    0.00287733,'AbsTol', 1e-8, 'Mthetatheta2 not behaving as expected')
end

function testMthetatheta2NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./2;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mthetatheta2NB(theta, phi, A, B, C, H, NF, a11, a13, a33),...
    0.00287646,'AbsTol', 1e-8, 'Mthetatheta2NB not behaving as expected')
end

function testQthetatheta3s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./2;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
Beta = 1e4;
I = 5000;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qthetatheta3(theta, phi, Beta, I, A, B, C, H, a11, a13, a33),...
    1.16639,'AbsTol', 1e-5, 'Qthetatheta3 not behaving as expected')
end

function testQthetatheta3NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./2;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
I = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qthetatheta3NB(theta, phi, I, A, B, C, H, a11, a13, a33),...
    2.33325,'AbsTol', 1e-5, 'Qthetatheta3NB not behaving as expected')
end

function testMthetatheta3s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./2;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
Beta = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mthetatheta3(theta, phi, Beta, A, B, C, H, NF, a11, a13, a33),...
    0.0140872,'AbsTol', 1e-7, 'Mthetatheta3 not behaving as expected')
end

function testMthetatheta3NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./2;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mthetatheta3NB(theta, phi, A, B, C, H, NF, a11, a13, a33),...
    0.014084,'AbsTol', 1e-6, 'Mthetatheta3NB not behaving as expected')
end

function testQthetatheta4s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
Beta = 1e4;
I = 5000;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qthetatheta4(theta, phi, Beta, I, A, B, C, H, a11, a13, a33),...
    0.218698,'AbsTol', 1e-6, 'Qthetatheta4 not behaving as expected')
end

function testQthetatheta4NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
I = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qthetatheta4NB(theta, phi, I, A, B, C, H, a11, a13, a33),...
    0.437484,'AbsTol', 1e-6, 'Qthetatheta4NB not behaving as expected')
end

function testMthetatheta4s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./2;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
Beta = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mthetatheta4(theta, phi, Beta, A, B, C, H, NF, a11, a13, a33),...
    0.00287733,'AbsTol', 1e-8, 'Mthetatheta4 not behaving as expected')
end

function testMthetatheta4NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./2;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mthetatheta4NB(theta, phi, A, B, C, H, NF, a11, a13, a33),...
    0.00287646,'AbsTol', 1e-8, 'Mthetatheta4NB not behaving as expected')
end

function testQthetaphi1s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
Beta = 1e4;
I = 500;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qthetaphi1(theta, phi, Beta, I, A, B, C, H, a11, a13, a33),...
    0.0857901,'AbsTol', 1e-7, 'Qthetaphi1 not behaving as expected')
end

function testQthetaphi1NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
I = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qthetaphi1NB(theta, phi, I, A, B, C, H, a11, a13, a33),...
    1.71615,'AbsTol', 1e-5, 'Qthetaphi1NB not behaving as expected')
end

function testMthetaphi1s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
Beta = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mthetaphi1(theta, phi, Beta, A, B, C, H, NF, a11, a13, a33),...
    0.0110329,'AbsTol', 1e-7, 'Mthetaphi1 not behaving as expected')
end

function testMthetaphi1NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mthetaphi1NB(theta, phi, A, B, C, H, NF, a11, a13, a33),...
    0.0110291,'AbsTol', 1e-7, 'Mthetaphi1NB not behaving as expected')
end

function testQthetaphi2s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
Beta = 1e4;
I = 500;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qthetaphi2(theta, phi, Beta, I, A, B, C, H, a11, a13, a33),...
    0.0428951,'AbsTol', 1e-7, 'Qthetaphi2 not behaving as expected')
end

function testQthetaphi2NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
I = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qthetaphi2NB(theta, phi, I, A, B, C, H, a11, a13, a33),...
    0.858073,'AbsTol', 1e-6, 'Qthetaphi2NB not behaving as expected')
end

function testMthetaphi2s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
Beta = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mthetaphi2(theta, phi, Beta, A, B, C, H, NF, a11, a13, a33),...
    0.00741719,'AbsTol', 1e-8, 'Mthetaphi2 not behaving as expected')
end

function testMthetaphi2NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mthetaphi2NB(theta, phi, A, B, C, H, NF, a11, a13, a33),...
    0.00741545,'AbsTol', 1e-8, 'Mthetaphi2NB not behaving as expected')
end

function testQthetaphi3s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
Beta = 1e4;
I = 500;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qthetaphi3(theta, phi, Beta, I, A, B, C, H, a11, a13, a33),...
    -0.0857901,'AbsTol', 1e-7, 'Qthetaphi3 not behaving as expected')
end

function testQthetaphi3NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
I = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qthetaphi3NB(theta, phi, I, A, B, C, H, a11, a13, a33),...
    -1.71615,'AbsTol', 1e-5, 'Qthetaphi3NB not behaving as expected')
end

function testMthetaphi3s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
Beta = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mthetaphi3(theta, phi, Beta, A, B, C, H, NF, a11, a13, a33),...
    0.0138158,'AbsTol', 1e-7, 'Mthetaphi3 not behaving as expected')
end

function testMthetaphi3NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mthetaphi3NB(theta, phi, A, B, C, H, NF, a11, a13, a33),...
    0.0138122,'AbsTol', 1e-7, 'Mthetaphi3NB not behaving as expected')
end

function testQthetaphi4s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
Beta = 1e4;
I = 500;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qthetaphi4(theta, phi, Beta, I, A, B, C, H, a11, a13, a33),...
    -0.0428951,'AbsTol', 1e-7, 'Qthetaphi4 not behaving as expected')
end

function testQthetaphi4NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
I = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qthetaphi4NB(theta, phi, I, A, B, C, H, a11, a13, a33),...
    -0.858073,'AbsTol', 1e-6, 'Qthetaphi4NB not behaving as expected')
end

function testMthetaphi4s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
Beta = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mthetaphi4(theta, phi, Beta, A, B, C, H, NF, a11, a13, a33),...
    0.00500718,'AbsTol', 1e-8, 'Mthetaphi4 not behaving as expected')
end

function testMthetaphi4NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mthetaphi4NB(theta, phi, A, B, C, H, NF, a11, a13, a33),...
    0.00500519,'AbsTol', 1e-8, 'Mthetaphi4NB not behaving as expected')
end

function testQphi1s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, ~, ~, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
Beta = 1e4;
I = 500;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qphi1(theta, phi, Beta, I, C, H, a11, a13, a33),...
    0.252401,'AbsTol', 1e-6, 'Qphi1 not behaving as expected')
end

function testQphi1NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, ~, ~, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
I = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qphi1NB(theta, phi, I, C, H, a11, a13, a33),...
    5.04903,'AbsTol', 1e-5, 'Qphi1NB not behaving as expected')
end

function testMphi1s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
Beta = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mphi1(theta, phi, Beta, A, B, C, H, NF, a11, a13, a33),...
    0.0119101,'AbsTol', 1e-7, 'Mphi1 not behaving as expected')
end

function testMphi1NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mphi1NB(theta, phi, A, B, C, H, NF, a11, a13, a33),...
    0.011906,'AbsTol', 1e-6, 'Mphi1NB not behaving as expected')
end

function testQphi2s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, ~, ~, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
Beta = 1e4;
I = 500;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qphi2(theta, phi, Beta, I, C, H, a11, a13, a33),...
    0.0841336,'AbsTol', 1e-7, 'Qphi2 not behaving as expected')
end

function testQphi2NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, ~, ~, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
I = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qphi2NB(theta, phi, I, C, H, a11, a13, a33),...
    1.68301,'AbsTol', 1e-5, 'Qphi2NB not behaving as expected')
end

function testMphi2s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
Beta = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mphi2(theta, phi, Beta, A, B, C, H, NF, a11, a13, a33),...
    0.0160138,'AbsTol', 1e-7, 'Mphi2 not behaving as expected')
end

function testMphi2NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mphi2NB(theta, phi, A, B, C, H, NF, a11, a13, a33),...
    0.0160101,'AbsTol', 1e-7, 'Mphi2NB not behaving as expected')
end

function testQphi3s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, ~, ~, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
Beta = 1e4;
I = 500;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qphi3(theta, phi, Beta, I, C, H, a11, a13, a33),...
    0.252401,'AbsTol', 1e-6, 'Qphi3 not behaving as expected')
end

function testQphi3NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, ~, ~, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
I = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qphi3NB(theta, phi, I, C, H, a11, a13, a33),...
    5.04903,'AbsTol', 1e-5, 'Qphi3NB not behaving as expected')
end

function testMphi3s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
Beta = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mphi3(theta, phi, Beta, A, B, C, H, NF, a11, a13, a33),...
    0.0149142,'AbsTol', 1e-7, 'Mphi3 not behaving as expected')
end

function testMphi3NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mphi3NB(theta, phi, A, B, C, H, NF, a11, a13, a33),...
    0.0149104,'AbsTol', 1e-7, 'Mphi3NB not behaving as expected')
end

function testQphi4s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, ~, ~, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
Beta = 1e4;
I = 500;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qphi4(theta, phi, Beta, I, C, H, a11, a13, a33),...
    0.0841336,'AbsTol', 1e-7, 'Qphi7 not behaving as expected')
end

function testQphi4NBs(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, ~, ~, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
I = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Qphi4NB(theta, phi, I, C, H, a11, a13, a33),...
    1.68301,'AbsTol', 1e-5, 'Qphi4NB not behaving as expected')
end

function testMphi4s(testCase)
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
theta = pi./4;
phi = pi./3;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
Beta = 1e4;
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mphi4(theta, phi, Beta, A, B, C, H, NF, a11, a13, a33),...
    0.0108106,'AbsTol', 1e-7, 'Mphi4 not behaving as expected')
verifyEqual(testCase, ...
    testCase.TestData.CRs.Mphi4NB(theta, phi, A, B, C, H, NF, a11, a13, a33),...
    0.0108063,'AbsTol', 1e-7, 'Mphi4NB not behaving as expected')
end

function testFIMs1(testCase)
theta = pi./3;
phi = pi./4;
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
time=1;
n_detectors = 4;
I = repmat(1000, n_detectors, 1);
Beta = repmat(1e4, n_detectors, 1);
FisherInfoM = testCase.TestData.CRs.FisherMatrixScatter(time, theta, phi, I, Beta,...
                A, B, C, H, NF, a11, a13, a33, n_detectors);
ExpcResult = [[13.730309080148039, 2.367497415301683e-16]; [2.367497415301683e-16, 1.588934160885147e+02]];
verifyEqual(testCase, FisherInfoM, ExpcResult, 'AbsTol', 1e-7, 'FIM (scattering) not behaving as expected')

theta = pi./2;
phi = pi./4;
NAObj = 0.7;
NACond = 1.3;
a11 = 0.3;
a13 = 0.1;
a33 = 0.6;
[~, ~, A, B, C, H] = testCase.TestData.CRs.InstrResp(NACond,NAObj,testCase.TestData.n0);
NF = testCase.TestData.CRs.NormFactor(theta, a11, a13, a33, A, B, H);
time=1;
n_detectors = 4;
I = repmat(1000, n_detectors, 1);
Beta = repmat(1e4, n_detectors, 1);
FisherInfoM = testCase.TestData.CRs.FisherMatrixScatter(time, theta, phi, I, Beta,...
                A, B, C, H, NF, a11, a13, a33, n_detectors);
ExpcResult = [[0, 0]; [0, 2.020738051470170E02]];
verifyEqual(testCase, FisherInfoM, ExpcResult, 'AbsTol', 1e-7, 'FIM (scattering) not behaving as expected')
end