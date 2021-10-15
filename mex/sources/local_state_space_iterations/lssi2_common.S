        ## Push a register to the stack and adjust CFI information accordingly
        ## Also increment a counter
        .macro push_cfi_reg reg
        push \reg
        .cfi_adjust_cfa_offset 8
        .cfi_rel_offset \reg, 0
        .set pushed_regs,pushed_regs+1
        .endm

        ## Pop a register from the stack and adjust CFI information accordingly
        .macro pop_cfi_reg reg
        pop \reg
        .cfi_adjust_cfa_offset -8
        .cfi_restore \reg
        .endm
