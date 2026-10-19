# NumericalSum_old.c → NumericalSum.c: what changed and why

`NumericalSum_old.c` (the first version) is unchanged. `NumericalSum.c` is the new version.

## Summary

| | `NumericalSum_old.c` (old) | `NumericalSum.c` (new) |
|---|---|---|
| **Input** | N and a from the command line. x and y typed in at prompts (or each element one by one). | Everything from the command line: `./numericalsum N a x y [expected]` |
| **`N = 1e8`** | `atoi("1e8")` reads only the `1`, so **N = 1**, silently | Read with `strtod`, so `1e8` = 100 000 000. Non-integer or negative N is rejected. |
| **Invalid input** (e.g. `abc`) | Silently becomes 0 | Error message + usage, exit code 1 |
| **Memory** | `malloc` not checked: crash if memory runs out | Checked. Prints the memory needed and a clear error if it fails. |
| **Output** | Prints **every** element (10⁸ lines ≈ 2.2 GB for N = 10⁸, ~13 s) | Prints the first 5 elements with 17 digits; the rest is only checked |
| **Point 3: check all elements** | ❌ Missing | ✅ Every `d[i]` compared with `a*x + y` |
| **Point 4: test = 7.4** | ❌ Missing | ✅ Optional `expected` argument, compared exactly (`==`) *and* with a tolerance |
| **Print precision** | `%f` (6 decimals): shows `7.400000` and hides the rounding error | `%.17g`: shows `7.3999999999999995` |
| **Timing** | None | Time of the `d = a*x + y` loop (`clock_gettime`) |
| **Exit code** | Always 0 | 0 = all OK, 1 = bad input / no memory, 2 = check failed |
| **Integer type for N** | `int` (max ≈ 2.1 × 10⁹) | `size_t` (the type C uses for sizes) |
| **Compiler warnings** (`-Wall -Wextra -O2`) | 8 warnings (unchecked `scanf`) | None |
| **Build** | Comment says `gcc NumericalSum.C` (capital `.C` = compiled as **C++**) | `Makefile`: just type `make` |

## The changes in more detail

### 1. Command-line input
The assignment asks for a, x and y "from command line". The old version asked for x and y interactively, which
makes the program impossible to run in a script or loop. In the new version, all values are arguments:

```bash
for N in 10 1e6 1e8; do ./numericalsum $N 3 0.1 7.1 7.4; done
```

### 2. Reading numbers safely
`atoi`/`atof` stop at the first character they don't understand and never report an error.
The new version uses `strtod` and checks that the *whole* string was a number:

```c
*out = strtod(s, &end);
return end != s && *end == '\0';   // something was read, and nothing is left over
```

### 3. Checking `malloc`
For N = 10⁸ the program needs 2.4 GB. If the system can't provide it, `malloc` returns `NULL`, and
the old version would crash writing to it. The new version checks:

```c
if (x == NULL || y == NULL || d == NULL) { fprintf(stderr, "Error: not enough memory ..."); return 1; }
```

### 4. Testing the result (points 3 and 4)
- **Check 1:** compares each `d[i]` with `a*x + y` computed once by the program. They are always
  identical: the same operations on the same numbers give the same rounding.
- **Check 2:** compares each `d[i]` with the value *you* expect (`7.4`):
  - `d[i] == 7.4` is **false**, because `3*0.1 + 7.1 = 7.3999999999999995` in floating point, while
    `7.4` is stored as `7.4000000000000004`;
  - `|d[i] − 7.4| ≤ 4·DBL_EPSILON·|7.4|` is **true**: they differ only by rounding (8.9 × 10⁻¹⁶).

  **Lesson:** never compare floating-point results with `==`. Use a tolerance.

### 5. Printing
The old version printed `%f` (6 decimals), so `7.3999999999999995` was displayed as `7.400000`. The rounding
error was invisible. The new version prints `%.17g`, which shows enough digits to identify the exact `double`.
It also prints only the first 5 elements, so N = 10⁸ runs in ~1 s instead of producing gigabytes
of text.
