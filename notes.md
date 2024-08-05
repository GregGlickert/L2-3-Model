currently the best performing model is the LTS2PNR4 block1. However the power of both beta and gamma is still fairly weak. You can see in the spike rate residuals that there is some beta and gamma, but i think it would need to be cleaned up and increased slightly still. The FSI2LTS and LTS2FSI syns could be helpful when it comes to this. Also during baseline the LTS is firing more than the FSI which i think may be a slight issue. 


UPDATE
LOVE the Pulse2PN block3 things are looking real solid now. Just need to calm down the FSI a bit in baseline and decrease PN firing also. Might be able to get away with just lowering base2PN. May need to also lower PN2FSI which is fine we have loads of gamma.

UPDATE 7/17
E rev for LTS/FSI2PN was set to -90mV this needs to be altered to -75mV. The Pulse2PN block3 still the best.

Best is now Thal2PN block2!

best is FSI2PNR2 block 7 looks awesome

slightly changed cell template and now model is firing a bit high for PN and FSI. Decreases PN2FSI from 3.1 to 2.9 and doing seedSweep on Thal2PN.


RESULT_PATH = "../Run-Storage/DONE_lower_FSI2PN_PulseSweepR2/block4"

RESULT_PATH = "../Run-Storage/reduce_FSI_FR_LTS2FSI/block3"

RESULT_PATH = "../Run-Storage/reduce_FSI_FR_LTS2FSI/block3"