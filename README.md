# BurstingBubble_VE_coated
 Bursting VE layer coated gas bubbles at a liquid-gas interface

## Use the stable version of the code:
- [burstingBubbles_viscoelastic_coated.c](burstingBubbles_viscoelastic_coated.c) where the coating is viscoelastic and wets both the bubble-bulk and the bulk-ambient surface (coating is fully wetting).
- [burstingBubbles_viscoelastic_coated_dewets.c](burstingBubbles_viscoelastic_coated_dewets.c) where the coating is viscoelastic but does not wet either the bubble-bulk or the bulk-ambient surface (coating is partially wetting). The coating will dewet from the bubble-bulk interface. Note that in this case, the bulk forms a precursor file. 

For details of this precursor film based model, see: [https://youtu.be/ozrnYe8u1HA?si=cH2FcO_kPw5o_6k7](https://youtu.be/ozrnYe8u1HA?si=cH2FcO_kPw5o_6k7)

## Known issues: 

- The fully elastic case [archived_codes/burstingBubbles_elastic_coated.c](archived_codes/burstingBubbles_elastic_coated.c) breaks down close to the singularity of focussing capillary waves. See [issue 15](https://github.com/VatsalSy/BurstingBubble_VE_coated/issues/15) for details. 


