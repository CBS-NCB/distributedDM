function trialSplits = generateTrialsSplits(dataBehaviorZSValidSet, curAnimalValidSet, validFrames, lateFrames, earlyFrames, presFrames, OLFrames)
  trialSplits = struct;
  %%% Wheel movements definitions
  % on zscored velocities, threshold of 2
  sigVel = sum(abs(dataBehaviorZSValidSet.wV(validFrames, :)) > 2);
  sigVelLate = sum(abs(dataBehaviorZSValidSet.wV(lateFrames, :)) > 2);
  sigVelOL = sum(abs(dataBehaviorZSValidSet.wV(OLFrames, :)) > 2);
  sigVelEarly = sum(abs(dataBehaviorZSValidSet.wV(earlyFrames, :)) > 2);
  sigVelP = sum((dataBehaviorZSValidSet.wV(validFrames, :)) > 2);
  sigVelPLate = sum((dataBehaviorZSValidSet.wV(lateFrames, :)) > 2);
  sigVelPEarly = sum((dataBehaviorZSValidSet.wV(earlyFrames, :)) > 2);
  sigVelN = sum((dataBehaviorZSValidSet.wV(validFrames, :)) < -2);
  sigVelNLate = sum((dataBehaviorZSValidSet.wV(lateFrames, :)) < -2);
  sigVelNEarly = sum((dataBehaviorZSValidSet.wV(earlyFrames, :)) < -2);
  
  trialSplits.all = 1:length(curAnimalValidSet.trialsSummary.choicedir);
  
  diffThres = [0 30 60 90];
  
  trialSplits.easy = find(abs(diff(abs(curAnimalValidSet.trialsSummary.ori),[],2)) <= diffThres(4) & abs(diff(abs(curAnimalValidSet.trialsSummary.ori),[],2)) > diffThres(3));
  trialSplits.med = find(abs(diff(abs(curAnimalValidSet.trialsSummary.ori),[],2)) <= diffThres(3) & abs(diff(abs(curAnimalValidSet.trialsSummary.ori),[],2)) > diffThres(2));
  trialSplits.hard = find(abs(diff(abs(curAnimalValidSet.trialsSummary.ori),[],2)) <= diffThres(2) & abs(diff(abs(curAnimalValidSet.trialsSummary.ori),[],2)) > diffThres(1));
  trialSplits.easyHalf = find(abs(diff(abs(curAnimalValidSet.trialsSummary.ori),[],2)) <= 90 & abs(diff(abs(curAnimalValidSet.trialsSummary.ori),[],2)) > 45);
  trialSplits.hardHalf = find(abs(diff(abs(curAnimalValidSet.trialsSummary.ori),[],2)) <= 45 & abs(diff(abs(curAnimalValidSet.trialsSummary.ori),[],2)) > 0);
  
  trialSplits.easiest = find(abs(diff(abs(curAnimalValidSet.trialsSummary.ori),[],2)) <= 90 & abs(diff(abs(curAnimalValidSet.trialsSummary.ori),[],2)) > 75);
  trialSplits.hardest = find(abs(diff(abs(curAnimalValidSet.trialsSummary.ori),[],2)) <= 15 & abs(diff(abs(curAnimalValidSet.trialsSummary.ori),[],2)) > 0);
  
  trialSplits.difficulty = abs(diff(abs(curAnimalValidSet.trialsSummary.ori),[],2));
  trialSplits.angleDiff = -diff(abs(curAnimalValidSet.trialsSummary.ori), [], 2);
  trialSplits.noWheel = find(sigVel < 5);
  trialSplits.wheel = find(sigVel >= 5);
  trialSplits.noWheelOL = find(sigVelOL < 5);
  trialSplits.WheelOL = find(sigVelOL < 5);

  trialSplits.noWheelLate = find(sigVelLate < 3); % 3 since there's less frames total
  trialSplits.wheelLate = find(sigVelLate >= 3);
  trialSplits.noWheelEarly = find(sigVelEarly < 5);
  trialSplits.wheelEarly = find(sigVelEarly >= 5);
  % Now with sign
  trialSplits.noWheelP = find(sigVelP < 5);
  trialSplits.wheelP = find(sigVelP >= 5);
  trialSplits.noWheelPLate = find(sigVelPLate < 3); % 3 since there's less frames total
  trialSplits.wheelPLate = find(sigVelPLate >= 3);
  trialSplits.noWheelPEarly = find(sigVelPEarly < 5);
  trialSplits.wheelPEarly = find(sigVelPEarly >= 5);
  
  trialSplits.noWheelN = find(sigVelN < 5);
  trialSplits.wheelN = find(sigVelN >= 5);
  trialSplits.noWheelNLate = find(sigVelNLate < 3); % 3 since there's less frames total
  trialSplits.wheelNLate = find(sigVelNLate >= 3);
  trialSplits.noWheelNEarly = find(sigVelNEarly < 5);
  trialSplits.wheelNEarly = find(sigVelNEarly >= 5);
  
  %%% Saccades definitions
  % Any saccade
  sigSacc = sum(abs(dataBehaviorZSValidSet.pS(validFrames, :)) > 0);
  sigSaccLate = sum(abs(dataBehaviorZSValidSet.pS(lateFrames, :)) > 0);
  sigSaccEarly = sum(abs(dataBehaviorZSValidSet.pS(earlyFrames, :)) > 0);
  sigSaccOL = sum(abs(dataBehaviorZSValidSet.pS(OLFrames, :)) > 0);
  
  trialSplits.noSacc = find(sigSacc == 0);
  trialSplits.sacc = find(sigSacc > 0);
  trialSplits.noSaccLate = find(sigSaccLate == 0);
  trialSplits.saccLate = find(sigSaccLate > 0);
  trialSplits.noSaccOL = find(sigSaccOL > 0);
  trialSplits.noSaccEarly = find(sigSaccEarly == 0);
  trialSplits.saccEarly = find(sigSaccEarly > 0);

  %%% pupX definitions
  sigPupX = mean(dataBehaviorZSValidSet.pX(presFrames, :));
  trialSplits.pupXcenter = find(abs(sigPupX) <= 1); % Around 1 std from center
  trialSplits.pupXnasal = find(sigPupX < -1);
  trialSplits.pupXtemporal = find(sigPupX > 1);
  
  %%% PupA definitions
  sigPupA = max(dataBehaviorZSValidSet.pA(OLFrames, :)-mean(dataBehaviorZSValidSet.pA(presFrames, :))); % Max OL - mean prestim
  trialSplits.pupAhigh = find(sigPupA >= prctile(sigPupA, 66));
  trialSplits.pupAlow = find(sigPupA <= prctile(sigPupA, 33));

  %%% History definitions (helpers)
  trialSplits.prevConsecutive = 1+find(diff(curAnimalValidSet.trialsSummary.trial) == 1);
  %trialSplits.prevPrevConsecutive = 2+find(curAnimalValidSet.trialsSummary.trial(3:end)-curAnimalValidSet.trialsSummary.trial(1:end-2) == 2);
  
  %%% Choice definitions
  trialSplits.left = find(curAnimalValidSet.trialsSummary.choicedir == -1);
  trialSplits.right = find(curAnimalValidSet.trialsSummary.choicedir == 1);
  trialSplits.timeout = find(curAnimalValidSet.trialsSummary.choicedir == 0);
  trialSplits.noTimeout = find(curAnimalValidSet.trialsSummary.choicedir ~= 0);
  trialSplits.prevLeft = intersect(find(curAnimalValidSet.trialsSummary.choicedir(1:end-1) == -1)+1, trialSplits.prevConsecutive);
  trialSplits.prevRight = intersect(find(curAnimalValidSet.trialsSummary.choicedir(1:end-1) == 1)+1, trialSplits.prevConsecutive);
  trialSplits.prevTimeout = intersect(find(curAnimalValidSet.trialsSummary.choicedir(1:end-1) == 0)+1, trialSplits.prevConsecutive);
  
  %%% Outocme definitions
  trialSplits.correct = find(curAnimalValidSet.trialsSummary.outcome == 1);
  trialSplits.incorrect = find(curAnimalValidSet.trialsSummary.outcome == -1);
  trialSplits.prevCorrect = intersect(find(curAnimalValidSet.trialsSummary.outcome(1:end-1) == 1)+1, trialSplits.prevConsecutive);
  trialSplits.prevIncorrect = intersect(find(curAnimalValidSet.trialsSummary.outcome(1:end-1) == -1)+1, trialSplits.prevConsecutive);
  
  %%% Strategies
  trialSplits.stay = union(intersect(trialSplits.prevLeft, trialSplits.left), intersect(trialSplits.prevRight, trialSplits.right));
  trialSplits.switch = union(intersect(trialSplits.prevLeft, trialSplits.right), intersect(trialSplits.prevRight, trialSplits.left));
end