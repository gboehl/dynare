/*
 * Local state-space iteration at order 2, without pruning.
 * Specialized kernel in x86-64 assembly using AVX-512F.
 * See the function prototype in lssi2.hh.
 *
 * WARNING: the number of particles must be a multiple of 8.
 *
 * TODO:
 * – Implement Windows ABI as an alternative
 * – Check whether requiring AVX-512VL (and thus replacing the 3 instructions
 *   which manipulate VSIB in zmm* by their equivalent that manipulate ymm*)
 *   improves performance. Actually, is there any CPU that implements AVX-512F
 *   but not AVX-512VL?
 * – Check whether using int64 inside VSIB (instead of int32) improves
 *   performance
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
        .arch .avx512f


### Some useful macros
#include "lssi2_common.S"


### Some pre-defined vector constants
        .section .rodata

        .align 8
double_0p5:
        .double 0.5

        .align 32
v8_int_0_1__7:
        .int 0, 1, 2, 3, 4, 5, 6, 7


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
        .global lssi2_avx512
lssi2_avx512:
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
        vpbroadcastd ymm15,xmm15                # ymm15 = [ int32: n (×8) ]
        vmovd xmm14,q
        vpbroadcastd ymm14,xmm14                # ymm14 = [ int32: q (×8) ]
        vmovd xmm13,m
        vpbroadcastd ymm13,xmm13                # ymm13 = [ int32: m (×8) ]

        vmovdqa ymm12,[rip+v8_int_0_1__7]
        vpmulld ymm12,ymm12,ymm15               # ymm12 = base VSIB for current particles in ŷ
        vmovdqa ymm11,[rip+v8_int_0_1__7]
        vpmulld ymm11,ymm11,ymm14               # ymm11 = base VSIB for current particles in ε
        vmovdqa ymm10,[rip+v8_int_0_1__7]
        vpmulld ymm10,ymm10,ymm13               # ymm10 = base VSIB for particles in y

        vpslld ymm15,ymm15,3                    # ymm15 = [ int32: 8n (×8) ]
        vpslld ymm14,ymm14,3                    # ymm14 = [ int32: 8q (×8) ]
        vpslld ymm13,ymm13,3                    # ymm13 = [ int32: 8m (×8) ]

        mov eax,1
        ## NB: We broadcast into zmm in order to avoid requiring AVX-512VL
        vpbroadcastd zmm8,eax                   # ymm8 = [ int32: 1 (×8) ]
        vbroadcastsd zmm16,[rip+double_0p5]     # zmm16 = [ double: 0.5 (×8) ]

        ## Enter the main loop
        xor ebx,ebx                             # ebx = particle counter

next_8_particles:
        xor ecx,ecx                             # ecx = variable counter
        vmovdqa ymm9,ymm10                      # ymm9 = VSIB for current particles in y

next_variable:
        ## Use ymm0 to store the 4 next-period particles for this specific variable
        ## Initialize ymm0 to the constant
        mov r8,[rsp+constant_offset]
        vbroadcastsd zmm0,[r8+rcx*8]

        ## Add ghx·ŷ
        mov ebp,n                               # ebp = number of remaining state variables
        vmovdqa ymm1,ymm12                      # ymm1 = VSIB for current particles in ŷ
        mov r8,[rsp+ghx_offset]
        lea r8,[r8+rcx*8]                       # r8 = pointer to ghxᵢⱼ
        .align 16                               # Recommendation p. 537 of Kusswurm (2018)
next_state:
        kxnorw k1,k1,k1
        vgatherdpd zmm2{k1},[yhat+ymm1*8]       # zmm2 = current particles for ŷⱼ
        vfmadd231pd zmm0,zmm2,[r8]{1to8}        # Add ghxᵢⱼ·ŷⱼ to current particles
        add r8,r11                              # Move to next column in ghx
        vpaddd ymm1,ymm1,ymm8                   # Update VSIB for next state (j=j+1)
        sub ebp,1
        jnz next_state

        ## Add ghu·ε
        mov ebp,q                               # ebp = number of remaining shocks
        vmovdqa ymm1,ymm11                      # ymm1 = VSIB for current particles in ε
        mov r8,[rsp+ghu_offset]
        lea r8,[r8+rcx*8]                       # r8 = pointer to ghuᵢⱼ
        .align 16
next_shock:
        kxnorw k1,k1,k1
        vgatherdpd zmm2{k1},[epsilon+ymm1*8]    # zmm2 = current particles for εⱼ
        vfmadd231pd zmm0,zmm2,[r8]{1to8}        # Add ghuᵢⱼ·εⱼ to current particles
        add r8,r11                              # Move to next column in ghu
        vpaddd ymm1,ymm1,ymm8                   # Update VSIB for next shock (j=j+1)
        sub ebp,1
        jnz next_shock

        ## Add ½ghxx·ŷ⊗ŷ
        xor ebp,ebp                             # Index of first state (j₁)
        mov r8,[rsp+ghxx_offset]
        lea r8,[r8+rcx*8]                       # r8 = pointer to ghxxᵢⱼ
        vmovdqa ymm1,ymm12                      # ymm1 = VSIB for current particles in ŷⱼ₁
        .align 16
next_state_state_1:
        mov r10d,ebp                            # Index of second state (j₂)
        vmovdqa ymm2,ymm1                       # ymm2 = VSIB for current particles in ŷⱼ₂
        .align 16
next_state_state_2:
        kxnorw k1,k1,k1
        vgatherdpd zmm4{k1},[yhat+ymm1*8]       # zmm4 = current particles for ŷⱼ₁
        kxnorw k2,k2,k2
        vgatherdpd zmm6{k2},[yhat+ymm2*8]       # zmm6 = current particles for ŷⱼ₂
        vmulpd zmm3,zmm4,zmm6                   # zmm3 = particles for ŷⱼ₁·ŷⱼ₂
        ## NB: We compare the zmm registers to avoid requiring AVX-512VL
        vpcmpeqd k3,zmm1,zmm2                   # k3[0:7]=1 if diagonal, 0 otherwise
        vmulpd zmm3{k3},zmm3,zmm16              # zmm3 = ½ŷⱼ₁·ŷⱼ₂ if diagonal, ŷⱼ₁·ŷⱼ₂ otherwise
        vfmadd231pd zmm0,zmm3,[r8]{1to8}        # Add (½?)ghxxᵢⱼ·ŷⱼ₁·ŷⱼ₂ to current particles
        vpaddd ymm2,ymm2,ymm8                   # Update VSIB for next second state (j₂=j₂+1)
        add r8,r11                              # Move to next column in ghxx
        add r10d,1
        cmp r10d,n
        jl next_state_state_2
        vpaddd ymm1,ymm1,ymm8                   # Update VSIB for next first state (j₁=j₁+1)
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
        vmovdqa ymm1,ymm11                      # ymm1 = VSIB for current particles in εⱼ₁
        .align 16
next_shock_shock_1:
        mov r10d,ebp                            # Index of second shock (j₂)
        vmovdqa ymm2,ymm1                       # ymm2 = VSIB for current particles in εⱼ₂
        .align 16
next_shock_shock_2:
        kxnorw k1,k1,k1
        vgatherdpd zmm4{k1},[epsilon+ymm1*8]    # zmm4 = current particles for εⱼ₁
        kxnorw k2,k2,k2
        vgatherdpd zmm6{k2},[epsilon+ymm2*8]    # zmm6 = current particles for εⱼ₂
        vmulpd zmm3,zmm4,zmm6                   # zmm3 = particles for εⱼ₁·εⱼ₂
        ## NB: We compare the zmm registers to avoid requiring AVX-512VL
        vpcmpeqd k3,zmm1,zmm2                   # k3[0:7]=1 if diagonal, 0 otherwise
        vmulpd zmm3{k3},zmm3,zmm16              # zmm3 = ½εⱼ₁εⱼ₂ if diagonal, εⱼ₁εⱼ₂ otherwise
        vfmadd231pd zmm0,zmm3,[r8]{1to8}        # Add ghuuᵢⱼ·εⱼ₁εⱼ₂ to current particles
        vpaddd ymm2,ymm2,ymm8                   # Update VSIB for next second shock (j₂=j₂+1)
        add r8,r11                              # Move to next column in ghuu
        add r10d,1
        cmp r10d,q
        jl next_shock_shock_2
        vpaddd ymm1,ymm1,ymm8                   # Update VSIB for next first shock (j₁=j₁+1)
        add ebp,1
        mov eax,ebp
        imul rax,r11
        add r8,rax                              # Jump several columns in ghuu
        cmp ebp,q
        jl next_shock_shock_1

        ## Add ghxu·ŷ⊗ε
        mov ebp,n                               # ebp = number of remaining states
        vmovdqa ymm1,ymm12                      # ymm1 = VSIB for current particles in ŷ
        mov r8,[rsp+ghxu_offset]
        lea r8,[r8+rcx*8]                       # r8 = pointer to ghxuᵢⱼ
        .align 16
next_state_2:
        mov eax,q                               # eax = number of remaining shocks
        vmovdqa ymm2,ymm11                      # ymm2 = VSIB for current particles in ε
        kxnorw k1,k1,k1
        vgatherdpd zmm4{k1},[yhat+ymm1*8]       # zmm4 = current particles for ŷⱼ₁
        .align 16
next_shock_2:
        kxnorw k2,k2,k2
        vgatherdpd zmm6{k2},[epsilon+ymm2*8]    # zmm6 = current particles for εⱼ₂
        vmulpd zmm6,zmm4,zmm6                   # zmm6 = particles for ŷⱼ₁εⱼ₂
        vfmadd231pd zmm0,zmm6,[r8]{1to8}        # Add ghxuᵢⱼ·ŷⱼ₁εⱼ₂ to current particles
        add r8,r11                              # Move to next column in ghxu
        vpaddd ymm2,ymm2,ymm8                   # Update VSIB for next shock (j₂=j₂+1)
        sub eax,1
        jnz next_shock_2
        vpaddd ymm1,ymm1,ymm8                   # Update VSIB for next state (j₁=j₁+1)
        sub ebp,1
        jnz next_state_2

        ## Save updated particles to memory
        kxnorw k1,k1,k1
        vscatterdpd [y+ymm9*8]{k1},zmm0

        ## Loop over variables
        vpaddd ymm9,ymm9,ymm8                   # Update VSIB for next variable in y
        add ecx,1
        cmp ecx,m
        jl next_variable

        ## Loop over particles
        add ebx,8
        vpaddd ymm12,ymm12,ymm15                # Update base VSIB for ŷ
        vpaddd ymm11,ymm11,ymm14                # Update base VSIB for ε
        vpaddd ymm10,ymm10,ymm13                # Update base VSIB for y
        cmp ebx,s
        jl next_8_particles

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
