#ifndef H_hugePages
#define H_hugePages

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#ifdef __linux__
    #include <sys/mman.h>
#endif

// Ask the kernel to back a large, not yet touched, allocation with transparent huge pages.
// The genome and suffix arrays are read at random positions, so with 4 KB pages most accesses
// also miss the TLB. This only changes how memory is mapped, never its contents.
// It is a hint: if THP is disabled ("never") or the call fails, nothing changes.
// Set STAR_NO_HUGEPAGES=1 in the environment to switch it off (for A/B timing).
inline void adviseHugePages(void *p, size_t len)
{
#if defined(__linux__) && defined(MADV_HUGEPAGE)
    const size_t minLen = (size_t)64 << 20; //only large arrays
    const uintptr_t hp = (uintptr_t)2 << 20; //2 MB, the PMD huge-page size on x86-64 and 4K-page aarch64
    if (p == NULL || len < minLen || std::getenv("STAR_NO_HUGEPAGES") != NULL)
        return;
    uintptr_t a = ((uintptr_t)p + hp - 1) & ~(hp - 1);
    uintptr_t b = ((uintptr_t)p + len) & ~(hp - 1);
    if (b > a)
        madvise((void *)a, b - a, MADV_HUGEPAGE);
#else
    (void)p; (void)len;
#endif
}

#endif
