#ifndef ALIGNCHECK_H
#define ALIGNCHECK_H

#include"../include/macro.h"

#include<stdlib.h>
#include<stdio.h>

inline void is_aligned(void const *p, size_t const byte_align, char *file, int const line)
	{
	REQUIRE((unsigned long) p % byte_align == 0, "error in alignment at file %s, line %d", file, line);
	}

#endif // ALIGNCHECK_H
