/*
 * Local state-space iteration at order 2, without pruning.
 * Specialized kernel in x86-64 assembly using AVX2 and FMA extensions.
 * See the function prototype in lssi2.hh.
 *
 * WARNING: the number of particles must be a multiple of 4.
 *
 * TODO:       
 * – Implement Windows ABI as an alternative
 */

/*
 * Copyright © 2022 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <https://www.gnu.org/licenses/>.
 */

#if defined(__x86_64__) && defined(__LP64__)

        .intel_syntax noprefix

        ## Enforce required CPU features, so that GNU as errors out if we use
        ## instructions not in that feature set
        .arch generic64
        .arch .avx2
        .arch .fma


### Some useful macros
#include "lssi2_common.S"

### Some pre-defined vector constants
        .section .rodata

        .align 8
double_1:
        .double 1.
double_0p5:
        .double 0.5

        .align 16
v4_int_0_1_2_3:
        .int 0, 1, 2, 3


### Layout of the current stack frame and of the beginning of the previous
### one (for arguments passed through the stack).
### This defines offset symbols through which variables and arguments can
### be conveniently referred.
        .struct 0
constant_offset:
        .space 8
ghu_offset:
        .space 8
ghx_offset:
        .space 8
saved_regs_offset:
        .space 6*8
return_address_offset:
        .space 8
        ## Beginning of the previous stack frame
ghxx_offset:
        .space 8
ghuu_offset:
        .space 8
ghxu_offset:
        .space 8
m_offset:
        .space 8
n_offset:
        .space 8
q_offset:
        .space 8
s_offset:
        .space 8


### Function body
        .text
        .global lssi2_avx2
lssi2_avx2:
        .cfi_startproc
        .set pushed_regs,0
        push_cfi_reg rbp
        push_cfi_reg rbx
        push_cfi_reg r12
        push_cfi_reg r13
        push_cfi_reg r14
        push_cfi_reg r15
        .if return_address_offset-saved_regs_offset != pushed_regs*8
        .error "Incorrect stack frame layout regarding saved registers"
        .endif

        ## Allocate space for variables in the current stack frame
        sub rsp,saved_regs_offset
        .cfi_adjust_cfa_offset saved_regs_offset

        ## Define more explicit names for registers that will remain associated
        ## to input arguments
        .set y,rdi
        .set yhat,rsi
        .set epsilon,rdx
        .set m,r12d
        .set n,r13d
        .set q,r14d
        .set s,r15d

        ## Initialize those registers that need it
        mov m,[rsp+m_offset]                    # number of states + observed vars
        mov n,[rsp+n_offset]                    # number of states
        mov q,[rsp+q_offset]                    # number of shocks
        mov s,[rsp+s_offset]                    # number of particles

        ## Conversely save some arguments to the stack
        mov [rsp+ghx_offset],rcx
        mov [rsp+ghu_offset],r8
        mov [rsp+constant_offset],r9

        ## Do nothing if s=0
        test s,s
        jz done

        ## Precompute 8m (used for moving pointers by one column inside y and gh* matrices)
        mov r11d,m
        shl r11,3                               # Use a 64-bit op, since m could theoretically be 2³¹-1

        ## Pre-compute some vectorized constants
        vmovd xmm15,n
        vpbroadcastd xmm15,xmm15                # xmm15 = [ int32: n (×4) ]
        vmovd xmm14,q
        vpbroadcastd xmm14,xmm14                # xmm14 = [ int32: q (×4) ]
        mov rax,0x8000000000000000
        vmovq xmm13,rax
        vpbroadcastq ymm13,xmm13                # ymm13 = [ int64: 0x8000000000000000 (×4) ]

        vmovdqa xmm12,[rip+v4_int_0_1_2_3]
        vpmulld xmm12,xmm12,xmm15               # xmm12 = base VSIB for current particles in ŷ
        vmovdqa xmm11,[rip+v4_int_0_1_2_3]
        vpmulld xmm11,xmm11,xmm14               # xmm11 = base VSIB for current particles in ε

        vpslld xmm15,xmm15,2                    # xmm15 = [ int32: 4n (×4) ]
        vpslld xmm14,xmm14,2                    # xmm14 = [ int32: 4q (×4) ]

        mov eax,1
        vmovd xmm10,eax
        vpbroadcastd xmm10,xmm10                # xmm10 = [ int32: 1 (×4) ]
        vbroadcastsd ymm9,[rip+double_1]        # ymm9 = [ double: 1. (×4) ]
        vbroadcastsd ymm8,[rip+double_0p5]      # ymm8 = [ double: 0.5 (×4) ]

        ## Enter the main loop
        xor ebx,ebx                             # ebx = particle counter

next_4_particles:
        xor ecx,ecx                             # ecx = variable counter

next_variable:
        ## Use ymm0 to store the 4 next-period particles for this specific variable
        ## Initialize ymm0 to the constant
        mov r8,[rsp+constant_offset]
        vbroadcastsd ymm0,[r8+rcx*8]

        ## Add ghx·ŷ
        mov ebp,n                               # ebp = number of remaining state variables
        vmovdqa xmm1,xmm12                      # xmm1 = VSIB for current particles in ŷ
        mov r8,[rsp+ghx_offset]
        lea r8,[r8+rcx*8]                       # r8 = pointer to ghxᵢⱼ
        .align 16                               # Recommendation p. 537 of Kusswurm (2018)
next_state:
        vmovdqa ymm4,ymm13                      # Mask for vgatherdpd (will be cleared by the instr.)
        vgatherdpd ymm2,[yhat+xmm1*8],ymm4      # ymm2 = current particles for ŷⱼ
        vbroadcastsd ymm3,[r8]                  # ymm3 = ghxᵢⱼ repeated 4 times
        vfmadd231pd ymm0,ymm2,ymm3              # Add ghxᵢⱼ·ŷⱼ to current particles
        add r8,r11                              # Move to next column in ghx
        vpaddd xmm1,xmm1,xmm10                  # Update VSIB for next state (j=j+1)
        sub ebp,1
        jnz next_state

        ## Add ghu·ε
        mov ebp,q                               # ebp = number of remaining shocks
        vmovdqa xmm1,xmm11                      # xmm1 = VSIB for current particles in ε
        mov r8,[rsp+ghu_offset]
        lea r8,[r8+rcx*8]                       # r8 = pointer to ghuᵢⱼ
        .align 16
next_shock:
        vmovdqa ymm4,ymm13
        vgatherdpd ymm2,[epsilon+xmm1*8],ymm4   # ymm2 = current particles for εⱼ
        vbroadcastsd ymm3,[r8]                  # ymm3 = ghuᵢⱼ repeated 4 times
        vfmadd231pd ymm0,ymm2,ymm3              # Add ghuᵢⱼ·εⱼ to current particles
        add r8,r11                              # Move to next column in ghu
        vpaddd xmm1,xmm1,xmm10                  # Update VSIB for next shock (j=j+1)
        sub ebp,1
        jnz next_shock

        ## Add ½ghxx·ŷ⊗ŷ
        xor ebp,ebp                             # Index of first state (j₁)
        mov r8,[rsp+ghxx_offset]
        lea r8,[r8+rcx*8]                       # r8 = pointer to ghxxᵢⱼ
        vmovdqa xmm1,xmm12                      # xmm1 = VSIB for current particles in ŷⱼ₁
        .align 16
next_state_state_1:
        mov r10d,ebp                            # Index of second state (j₂)
        vmovdqa xmm2,xmm1                       # xmm2 = VSIB for current particles in ŷⱼ₂
        .align 16
next_state_state_2:
        vmovdqa ymm3,ymm13
        vgatherdpd ymm4,[yhat+xmm1*8],ymm3      # ymm4 = current particles for ŷⱼ₁
        vmovdqa ymm5,ymm13
        vgatherdpd ymm6,[yhat+xmm2*8],ymm5      # ymm6 = current particles for ŷⱼ₂
        vmulpd ymm3,ymm4,ymm6                   # ymm3 = particles for ŷⱼ₁·ŷⱼ₂
        vpcmpeqd xmm7,xmm1,xmm2
        vinserti128 ymm7,ymm7,xmm7,1            # ymm7 = all 1s if diagonal element, all 0s otherwise
        ## There is a spurious dependency here due to ymm4 being reused (no more free registers)
        ## Freeing ymm8 by using a memory read in the following blend does not improve things
        vblendvpd ymm4,ymm9,ymm8,ymm7           # ymm4 = ½ if diagonal element, 1 otherwise
        vmulpd ymm3,ymm4,ymm3                   # ymm3 = ½ŷⱼ₁·ŷⱼ₂ if diagonal, ŷⱼ₁·ŷⱼ₂ otherwise
        vbroadcastsd ymm7,[r8]                  # ymm7 = ghxxᵢⱼ repeated 4 times
        vfmadd231pd ymm0,ymm7,ymm3              # Add (½?)ghxxᵢⱼ·ŷⱼ₁·ŷⱼ₂ to current particles
        vpaddd xmm2,xmm2,xmm10                  # Update VSIB for next second state (j₂=j₂+1)
        add r8,r11                              # Move to next column in ghxx
        add r10d,1
        cmp r10d,n
        jl next_state_state_2
        vpaddd xmm1,xmm1,xmm10                  # Update VSIB for next first state (j₁=j₁+1)
        add ebp,1
        mov eax,ebp
        imul rax,r11
        add r8,rax                              # Jump several columns in ghxx
        cmp ebp,n
        jl next_state_state_1

        ## Add ½ghuu·ε⊗ε
        xor ebp,ebp                             # Index of first shock (j₁)
        mov r8,[rsp+ghuu_offset]
        lea r8,[r8+rcx*8]                       # r8 = pointer to ghuuᵢⱼ
        vmovdqa xmm1,xmm11                      # xmm1 = VSIB for current particles in εⱼ₁
        .align 16
next_shock_shock_1:
        mov r10d,ebp                            # Index of second shock (j₂)
        vmovdqa xmm2,xmm1                       # xmm2 = VSIB for current particles in εⱼ₂
        .align 16
next_shock_shock_2:
        vmovdqa ymm3,ymm13
        vgatherdpd ymm4,[epsilon+xmm1*8],ymm3   # ymm4 = current particles for εⱼ₁
        vmovdqa ymm5,ymm13
        vgatherdpd ymm6,[epsilon+xmm2*8],ymm5   # ymm6 = current particles for εⱼ₂
        vmulpd ymm3,ymm4,ymm6                   # ymm3 = particles for εⱼ₁·εⱼ₂
        vpcmpeqd xmm7,xmm1,xmm2
        vinserti128 ymm7,ymm7,xmm7,1            # ymm7 = all 1s if diagonal element, all 0s otherwise
        ## Same remark as previous section regarding ymm4
        vblendvpd ymm4,ymm9,ymm8,ymm7           # ymm4 = ½ if diagonal element, 1 otherwise
        vmulpd ymm3,ymm4,ymm3                   # ymm3 = ½εⱼ₁εⱼ₂ if diagonal, εⱼ₁εⱼ₂ otherwise
        vbroadcastsd ymm7,[r8]                  # ymm7 = ghuuᵢⱼ repeated 4 times
        vfmadd231pd ymm0,ymm7,ymm3              # Add ghuuᵢⱼ·εⱼ₁εⱼ₂ to current particles
        vpaddd xmm2,xmm2,xmm10                  # Update VSIB for next second shock (j₂=j₂+1)
        add r8,r11                              # Move to next column in ghuu
        add r10d,1
        cmp r10d,q
        jl next_shock_shock_2
        vpaddd xmm1,xmm1,xmm10                  # Update VSIB for next first shock (j₁=j₁+1)
        add ebp,1
        mov eax,ebp
        imul rax,r11
        add r8,rax                              # Jump several columns in ghuu
        cmp ebp,q
        jl next_shock_shock_1

        ## Add ghxu·ŷ⊗ε
        mov ebp,n                               # ebp = number of remaining states
        vmovdqa xmm1,xmm12                      # xmm1 = VSIB for current particles in ŷ
        mov r8,[rsp+ghxu_offset]
        lea r8,[r8+rcx*8]                       # r8 = pointer to ghxuᵢⱼ
        .align 16
next_state_2:
        mov eax,q                               # eax = number of remaining shocks
        vmovdqa xmm2,xmm11                      # xmm2 = VSIB for current particles in ε
        vmovdqa ymm3,ymm13
        vgatherdpd ymm4,[yhat+xmm1*8],ymm3      # ymm4 = current particles for ŷⱼ₁
        .align 16
next_shock_2:
        vmovdqa ymm5,ymm13
        vgatherdpd ymm6,[epsilon+xmm2*8],ymm5   # ymm6 = current particles for εⱼ₂
        vbroadcastsd ymm7,[r8]                  # ymm7 = ghxuᵢⱼ repeated 4 times
        vmulpd ymm6,ymm4,ymm6                   # ymm6 = particles for ŷⱼ₁εⱼ₂
        vfmadd231pd ymm0,ymm7,ymm6              # Add ghxuᵢⱼ·ŷⱼ₁εⱼ₂ to current particles
        add r8,r11                              # Move to next column in ghxu
        vpaddd xmm2,xmm2,xmm10                  # Update VSIB for next shock (j₂=j₂+1)
        sub eax,1
        jnz next_shock_2
        vpaddd xmm1,xmm1,xmm10                  # Update VSIB for next state (j₁=j₁+1)
        sub ebp,1
        jnz next_state_2

        ## Save updated particles to memory
        mov eax,ebx
        imul eax,m
        add eax,ecx
        lea r8,[y+rax*8]
        vextractf128 xmm2,ymm0,1
        vmovsd [r8],xmm0
        add r8,r11
        vpsrldq xmm1,xmm0,8
        vmovsd [r8],xmm1
        add r8,r11
        vmovsd [r8],xmm2
        add r8,r11
        vpsrldq xmm3,xmm2,8
        vmovsd [r8],xmm3

        ## Loop over variables
        add ecx,1
        cmp ecx,m
        jl next_variable

        ## Loop over particles
        add ebx,4
        vpaddd xmm12,xmm12,xmm15                # Update base VSIB for ŷ
        vpaddd xmm11,xmm11,xmm14                # Update base VSIB for ε
        cmp ebx,s
        jl next_4_particles

done:
        ## Cleanup
        add rsp,saved_regs_offset
        .cfi_adjust_cfa_offset -saved_regs_offset
        pop_cfi_reg r15
        pop_cfi_reg r14
        pop_cfi_reg r13
        pop_cfi_reg r12
        pop_cfi_reg rbx
        pop_cfi_reg rbp
        vzeroupper
        ret
        .cfi_endproc

#endif
