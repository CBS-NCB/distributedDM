function [curSplit, curSplitLabels, curSplitTitle] = chooseSplit(trialSplits, selectedSplit)
  switch selectedSplit
    case 'all-LR'
      curSplit = {trialSplits.left, ...
                 trialSplits.right};
      curSplitLabels = {'left', 'right'};
      curSplitTitle = 'All - L vs R';
    case 'all-EH'
      curSplit = {trialSplits.easy, ...
                 trialSplits.hard};
      curSplitLabels = {'easy', 'hard'};
      curSplitTitle = 'All - E vs H';
    case 'all-EHest'
      curSplit = {trialSplits.easiest, ...
                 trialSplits.hardest};
      curSplitLabels = {'easiest', 'hardest'};
      curSplitTitle = 'All - E vs H est';
    case 'L-EH'
      curSplit = {trialSplits.easy, ...
                 trialSplits.hard};
      curSplit = cellfun(@(x)intersect(x, trialSplits.left), curSplit, 'UniformOutput', false);
      curSplitLabels = {'easy', 'hard'};
      curSplitTitle = 'Left - E vs H';
    case 'CL-EHest'
      curSplit = {trialSplits.easiest, ...
                 trialSplits.hardest};
      curSplit = cellfun(@(x)intersect(x, trialSplits.left), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.correct), curSplit, 'UniformOutput', false);
      curSplitLabels = {'easiest', 'hardest'};
      curSplitTitle = 'C Left - Eest vs Hest';
    case 'CR-EHest'
      curSplit = {trialSplits.easiest, ...
                 trialSplits.hardest};
      curSplit = cellfun(@(x)intersect(x, trialSplits.right), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.correct), curSplit, 'UniformOutput', false);
      curSplitLabels = {'easiest', 'hardest'};
      curSplitTitle = 'C Right  - Eest vs Hest';
    case 'R-EH'
      curSplit = {trialSplits.easy, ...
                 trialSplits.hard};
      curSplit = cellfun(@(x)intersect(x, trialSplits.right), curSplit, 'UniformOutput', false);
      curSplitLabels = {'easy', 'hard'};
      curSplitTitle = 'Right - E vs H';
    case 'all-EHhalf'
      curSplit = {trialSplits.easyHalf, ...
                 trialSplits.hardHalf};
      curSplitLabels = {'easy half', 'hard half'};
      curSplitTitle = 'All - E vs H half';
    case 'Ceasy-LR'
      curSplit = {trialSplits.left, trialSplits.right};
      curSplit = cellfun(@(x)intersect(x, trialSplits.correct), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.easy), curSplit, 'UniformOutput', false);
      curSplitLabels = {'left Ceasy', 'right Ceasy'};
      curSplitTitle = 'Ceasy - L vs R';
    case 'Cmed-LR'
      curSplit = {trialSplits.left, trialSplits.right};
      curSplit = cellfun(@(x)intersect(x, trialSplits.correct), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.med), curSplit, 'UniformOutput', false);
      curSplitLabels = {'left Cmed', 'right Cmed'};
      curSplitTitle = 'Cmed - L vs R';
    case 'Chard-LR'
      curSplit = {trialSplits.left, trialSplits.right};
      curSplit = cellfun(@(x)intersect(x, trialSplits.correct), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.hard), curSplit, 'UniformOutput', false);
      curSplitLabels = {'left Chard', 'right Chard'};
      curSplitTitle = 'Chard - L vs R';
    case 'easy-LR'
      curSplit = {trialSplits.left, trialSplits.right};
      curSplit = cellfun(@(x)intersect(x, trialSplits.easy), curSplit, 'UniformOutput', false);
      curSplitLabels = {'left easy', 'right easy'};
      curSplitTitle = 'easy - L vs R';
    case 'med-LR'
      curSplit = {trialSplits.left, trialSplits.right};
      curSplit = cellfun(@(x)intersect(x, trialSplits.med), curSplit, 'UniformOutput', false);
      curSplitLabels = {'left med', 'right med'};
      curSplitTitle = 'med - L vs R';
    case 'hard-LR'
      curSplit = {trialSplits.left, trialSplits.right};
      curSplit = cellfun(@(x)intersect(x, trialSplits.hard), curSplit, 'UniformOutput', false);
      curSplitLabels = {'left hard', 'right hard'};
      curSplitTitle = 'hard - L vs R';
    case 'rand-LR'
      l1 = sort(union(trialSplits.left, trialSplits.right));
      l2 = l1(randperm(length(l1)));
      curSplit = {l2(1:length(trialSplits.left)), ...
                 l2((length(trialSplits.left)+1):end)};
      curSplitLabels = {'left', 'right'};
      curSplitTitle = 'Rand - L vs R';
    case 'C-LR'
      curSplit = {intersect(trialSplits.left, trialSplits.correct)', ...
                 intersect(trialSplits.right, trialSplits.correct)'};
      curSplitLabels = {'left correct', 'right correct'};
      curSplitTitle = 'Correct - L vs R';
    case 'CnoWnoS-LR'
      curSplit = {trialSplits.left, trialSplits.right};
      curSplit = cellfun(@(x)intersect(x, trialSplits.correct), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.noWheel), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.noSacc), curSplit, 'UniformOutput', false);
      curSplitLabels = {'left CnoWnoS', 'right CnoWnoS'};
      curSplitTitle = 'CnoWnoS - L vs R';
    case 'CEnoWnoS-LR'
      curSplit = {trialSplits.left, trialSplits.right};
      curSplit = cellfun(@(x)intersect(x, trialSplits.correct), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.easy), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.noWheel), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.noSacc), curSplit, 'UniformOutput', false);
      curSplitLabels = {'left CEnoWnoS', 'right CEnoWnoS'};
      curSplitTitle = 'CEnoWnoS - L vs R';
    case 'CnoWOLnoSOL-LR'
      curSplit = {trialSplits.left, trialSplits.right};
      curSplit = cellfun(@(x)intersect(x, trialSplits.correct), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.noWheelOL), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.noSaccOL), curSplit, 'UniformOutput', false);
      curSplitLabels = {'left CnoWOLnoSOL', 'right CnoWOLnoSOL'};
      curSplitTitle = 'CnoWOLnoSOL - L vs R';
    case 'noWOLnoSOL-LR'
      curSplit = {trialSplits.left, trialSplits.right};
      curSplit = cellfun(@(x)intersect(x, trialSplits.noWheelOL), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.noSaccOL), curSplit, 'UniformOutput', false);
      curSplitLabels = {'left noWOLnoSOL', 'right noWOLnoSOL'};
      curSplitTitle = 'noWOLnoSOL - L vs R';
    case 'CnoWOL-LR'
      curSplit = {trialSplits.left, trialSplits.right};
      curSplit = cellfun(@(x)intersect(x, trialSplits.correct), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.noWheelOL), curSplit, 'UniformOutput', false);
      curSplitLabels = {'left CnoWOL', 'right CnoWOL'};
      curSplitTitle = 'CnoWOL - L vs R';
    case 'CnoWnoSpX0-LR'
      curSplit = {trialSplits.left, trialSplits.right};
      curSplit = cellfun(@(x)intersect(x, trialSplits.correct), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.noWheel), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.noSacc), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.pupXcenter), curSplit, 'UniformOutput', false);
      curSplitLabels = {'left CnoWnoSpX0', 'right CnoWnoSpX0'};
      curSplitTitle = 'CnoWnoSpX0 - L vs R';
      
    case 'I-LR'
      curSplit = {intersect(trialSplits.left, trialSplits.incorrect), ...
                 intersect(trialSplits.right, trialSplits.incorrect)};
      curSplitLabels = {'left incorrect', 'right incorrect'};
      curSplitTitle = 'Inorrect - L vs R';
    case 'L-IC'
      curSplit = {intersect(trialSplits.left, trialSplits.incorrect), ...
                 intersect(trialSplits.left, trialSplits.correct)};
      curSplitLabels = {'left incorrect', 'left correct'};
      curSplitTitle = 'Left - I vs C';
    case 'R-IC'
      curSplit = {intersect(trialSplits.right, trialSplits.incorrect), ...
                 intersect(trialSplits.right, trialSplits.correct)};
      curSplitLabels = {'right incorrect', 'right correct'};
      curSplitTitle = 'Right - I vs C';
    case 'C-prevIC'
      curSplit = {trialSplits.prevIncorrect, trialSplits.prevCorrect};
      curSplit = cellfun(@(x)intersect(x, trialSplits.correct), curSplit, 'UniformOutput', false);
      curSplitLabels = {'C - prevI', 'C - prevC'};
      curSplitTitle = 'Correct - prevI vs prevC';
    case 'all-prevIC'
      curSplit = {trialSplits.prevIncorrect, trialSplits.prevCorrect};
      curSplitLabels = {'prevI', 'prevC'};
      curSplitTitle = 'All - prevI vs prevC';
    case 'all-pupALH'
      curSplit = {trialSplits.pupAlow, trialSplits.pupAhigh};
      curSplitLabels = {'pupA low', 'pupA high'};
      curSplitTitle = 'All - pupA low vs pupA high';
    case 'all-lateWheel'
      curSplit = {trialSplits.noWheel, trialSplits.wheelLate};
      curSplitLabels = {'no wheel', 'wheel late'};
      curSplitTitle = 'All - late wheel vs no wheel';
    case 'CnoSaccNoEarlyWheel-lateWheel'
      curSplit = {trialSplits.noWheel, trialSplits.wheelLate};
      curSplit = cellfun(@(x)intersect(x, trialSplits.noWheelEarly), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.noSacc), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.correct), curSplit, 'UniformOutput', false);
      curSplitLabels = {'no wheel', 'wheel late'};
      curSplitTitle = 'C no Sacc - only late wheel vs no wheel';
    case 'CnoWnoSaccEarly-lateSacc'
      curSplit = {trialSplits.noSacc, trialSplits.saccLate};
      curSplit = cellfun(@(x)intersect(x, trialSplits.noWheel), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.noSaccEarly), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.correct), curSplit, 'UniformOutput', false);
      curSplitLabels = {'no sacc', 'sacc late'};
      curSplitTitle = 'C no W - only late sacc vs no sacc';
    case 'noSaccNoEarlyWheel-lateWheel'
      curSplit = {trialSplits.noWheel, trialSplits.wheelLate};
      curSplit = cellfun(@(x)intersect(x, trialSplits.noWheelEarly), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.noSacc), curSplit, 'UniformOutput', false);
      curSplitLabels = {'no wheel', 'wheel late'};
      curSplitTitle = 'no Sacc - only late wheel vs no wheel';
    case 'noWnoSaccEarly-lateSacc'
      curSplit = {trialSplits.noSacc, trialSplits.saccLate};
      curSplit = cellfun(@(x)intersect(x, trialSplits.noWheel), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.noSaccEarly), curSplit, 'UniformOutput', false);
      curSplitLabels = {'no sacc', 'sacc late'};
      curSplitTitle = 'no W - only late sacc vs no sacc';
    case 'CnoSaccNoEarylWheelLateWheel-LR'
      curSplit = {trialSplits.left, trialSplits.right};
      curSplit = cellfun(@(x)intersect(x, trialSplits.wheelLate), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.noWheelEarly), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.noSacc), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.correct), curSplit, 'UniformOutput', false);
      curSplitLabels = {'left', 'right'};
      curSplitTitle = 'C no Sacc only late wheel - left vs right';
    case 'CnoSaccNoEarylWheelLateWheel-wNwP'
      curSplit = {trialSplits.wheelNLate, trialSplits.wheelPLate};
      curSplit = cellfun(@(x)intersect(x, trialSplits.wheelLate), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.noWheelEarly), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.noSacc), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.correct), curSplit, 'UniformOutput', false);
      curSplitLabels = {'WN', 'WP'};
      curSplitTitle = 'C no Sacc only late wheel - WP vs WN';
    case 'noSaccNoEarylWheel-lateWheel'
      curSplit = {trialSplits.noWheel, trialSplits.wheelLate};
      curSplit = cellfun(@(x)intersect(x, trialSplits.noWheelEarly), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.noSacc), curSplit, 'UniformOutput', false);
      curSplitLabels = {'no wheel', 'wheel late'};
      curSplitTitle = 'no Sacc - only late wheel vs no wheel';
    case 'all-wheelPN'
      curSplit = {trialSplits.wheelP, trialSplits.wheelN};
      curSplitLabels = {'wheel P', 'wheel N'};
      curSplitTitle = 'All - wheel P vs wheel N';
    case 'all-lateWheelPN'
      curSplit = {trialSplits.wheelPLate, trialSplits.wheelNLate};
      curSplitLabels = {'wheel late P', 'wheel late N'};
      curSplitTitle = 'All - late wheel P vs late wheel N';
    case 'CpupAH-LR'
      curSplit = {trialSplits.left, trialSplits.right};
      curSplit = cellfun(@(x)intersect(x, trialSplits.correct), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.pupAhigh), curSplit, 'UniformOutput', false);
      curSplitLabels = {'left CpupAH', 'right CpupAH'};
      curSplitTitle = 'CpupAH - L vs R';
    case 'CpupAL-LR'
      curSplit = {trialSplits.left, trialSplits.right};
      curSplit = cellfun(@(x)intersect(x, trialSplits.correct), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.pupAlow), curSplit, 'UniformOutput', false);
      curSplitLabels = {'left CpupAL', 'right CpupAL'};
      curSplitTitle = 'CpupAL - L vs R';
    case 'CprevC-LR'
      curSplit = {trialSplits.left, trialSplits.right};
      curSplit = cellfun(@(x)intersect(x, trialSplits.correct), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.prevCorrect), curSplit, 'UniformOutput', false);
      curSplitLabels = {'left CprevC', 'right CprevC'};
      curSplitTitle = 'CprevC - L vs R';
    case 'CprevI-LR'
      curSplit = {trialSplits.left, trialSplits.right};
      curSplit = cellfun(@(x)intersect(x, trialSplits.correct), curSplit, 'UniformOutput', false);
      curSplit = cellfun(@(x)intersect(x, trialSplits.prevIncorrect), curSplit, 'UniformOutput', false);
      curSplitLabels = {'left CprevI', 'right CprevI'};
      curSplitTitle = 'CprevI - L vs R';
    case 'pupAH-LR'
      curSplit = {trialSplits.left, trialSplits.right};
      curSplit = cellfun(@(x)intersect(x, trialSplits.pupAhigh), curSplit, 'UniformOutput', false);
      curSplitLabels = {'left pupAH', 'right pupAH'};
      curSplitTitle = 'pupAH - L vs R';
    case 'pupAL-LR'
      curSplit = {trialSplits.left, trialSplits.right};
      curSplit = cellfun(@(x)intersect(x, trialSplits.pupAlow), curSplit, 'UniformOutput', false);
      curSplitLabels = {'left pupAL', 'right pupAL'};
      curSplitTitle = 'pupAL - L vs R';
    case 'prevC-LR'
      curSplit = {trialSplits.left, trialSplits.right};
      curSplit = cellfun(@(x)intersect(x, trialSplits.prevCorrect), curSplit, 'UniformOutput', false);
      curSplitLabels = {'left prevC', 'right prevC'};
      curSplitTitle = 'prevC - L vs R';
    case 'prevI-LR'
      curSplit = {trialSplits.left, trialSplits.right};
      curSplit = cellfun(@(x)intersect(x, trialSplits.prevIncorrect), curSplit, 'UniformOutput', false);
      curSplitLabels = {'left prevI', 'right prevI'};
      curSplitTitle = 'prevI - L vs R';
    case 'SwSt'
      curSplit = {trialSplits.switch, trialSplits.stay};
      curSplitLabels = {'switch', 'stay'};
      curSplitTitle = 'All - switch vs stay';
    case 'prevC-SwSt'
      curSplit = {intersect(trialSplits.switch, trialSplits.prevCorrect), intersect(trialSplits.stay, trialSplits.prevCorrect)};
      curSplitLabels = {'win switch', 'win stay'};
      curSplitTitle = 'prevC (win) - switch vs stay';
    case 'prevCpupAH-SwSt'
      curSplit = {intersect(intersect(trialSplits.switch, trialSplits.prevCorrect), trialSplits.pupAhigh), 
                  intersect(intersect(trialSplits.stay, trialSplits.prevCorrect), trialSplits.pupAhigh)};
      curSplitLabels = {'win switch', 'win stay'};
      curSplitTitle = 'prevC (win) pupAhigh - switch vs stay';
    case 'prevI-SwSt'
      curSplit = {intersect(trialSplits.switch, trialSplits.prevIncorrect), intersect(trialSplits.stay, trialSplits.prevIncorrect)};
      curSplitLabels = {'lose switch', 'lose stay'};
      curSplitTitle = 'prevI (lose) - switch vs stay';
    otherwise
      curSplit = [];
      curSplitLabels = {};
      curSplitTitle = 'invalid split';
  end
end
