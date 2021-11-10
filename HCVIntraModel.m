function dx = HCVIntraModel(t,x)
dx = zeros(17,1);

kTC = 1; ktrans = 180; kcleavage = 9;

kHFC = 0.0008; kRp5B = 0.1; kRip = 0.6; kinit = 1.12; kRids = 10; krepl = 1.12; 

koutRp = 0.307; kassembly = 1.2e-7;

kdegRp = 0.26;  
if (t<=21)
    kdegS = 0.61;
else
    kdegS = 0.1;
end

kdegNS = 0.11; kdegVMS = 0.001;

ks = kcleavage;

dx(1) = -kTC*x(2)*x(1) + ktrans*x(3) - kRp5B*x(1)*x(8) + koutRp*x(16) - kdegRp*x(1);
dx(2) = -kTC*x(2)*x(1) + ktrans*x(3);
dx(3) = kTC*x(2)*x(1) - ktrans*x(3);
dx(4) = ktrans*x(3) - kcleavage*x(4);
dx(5) = ks*x(4) - kdegS*x(5) - 180*kassembly*x(16)*x(5);
dx(6) = kcleavage*x(4) - kdegNS*x(6);
dx(7) = kcleavage*x(4) - kdegNS*x(7) - kHFC*x(7)*x(9);
dx(8) = kcleavage*x(4) - kdegNS*x(8) - kRp5B*x(1)*x(8);
dx(9) = -kHFC*x(7)*x(9) + kinit*x(12);
dx(10) = kHFC*x(7)*x(9) - kRip*x(11)*x(10);
dx(11) = kRp5B*x(1)*x(8) - kRip*x(11)*x(10);
dx(12) = kRip*x(11)*x(10) - kinit*x(12) - kdegVMS*x(12);
dx(13) = kinit*x(12) - kRids*x(13)*x(14) - kdegVMS*x(13) + krepl*x(15);
dx(14) = kinit*x(12) - kRids*x(13)*x(14) - kdegVMS*x(14) + krepl*x(15);
dx(15) = kRids*x(13)*x(14) - kdegVMS*x(15) - krepl*x(15);
dx(16) = krepl*x(15) - kassembly*x(16)*x(5) - koutRp*x(16) - kdegVMS*x(16);
dx(17) = kassembly*x(16)*x(5);
end
