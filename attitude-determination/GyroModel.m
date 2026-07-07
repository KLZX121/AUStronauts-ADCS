function [omegaMeas, d] = GyroModel(d, omegaTrue, dt)

d = RateGyro(d, omegaTrue, dt);

omegaMeas = zeros(3,1);

for i = 1:3
    omegaMeas(i) = d(i).measurement * d(i).resolution;
end

end