#include "cpucycles.h"

/**
 * Returns the current value of the CPU's time-stamp counter.
 * Uses inline assembly to access the RDTSC (Read Time-Stamp Counter) instruction,
 * which reads the time-stamp counter of the CPU into the EDX:EAX registers.
 *
 * @return: The current value of the time-stamp counter.
 */
long long cpucycles()
{
    unsigned long long result;  // Variable to store the result of the RDTSC instruction.

    // Inline assembly to execute the RDTSC instruction.
    // ".byte 15;.byte 49;" is an older way to represent the RDTSC instruction.
    // "shlq $32,%%rdx" shifts the value in the rdx register left by 32 bits.
    // "orq %%rdx,%%rax" performs a bitwise OR between rax and rdx, storing the result in rax.
    // The final 64-bit timestamp counter value is composed of the higher 32 bits from rdx and the lower 32 bits from rax.
    asm volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax"
            : "=a" (result)  // Output operand: The result is stored in the 'result' variable using the rax register.
            :                // No input operands.
            :  "%rdx"        // Clobbered register: rdx is modified by the assembly code.
            );

    return result;  // Return the 64-bit timestamp counter value.
}
