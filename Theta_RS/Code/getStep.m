function Ladder = getStep(T, dt, noSteps, ton)

t = 0:dt:T;

stepLength = T/noSteps;

LastStep = floor(t./stepLength)';

StepT = mod(t, stepLength);

Step = ((StepT/ton).*(StepT<=ton)+(ton<StepT))';

Ladder = Step + LastStep;