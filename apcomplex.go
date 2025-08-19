// Package apcomplex provides arbitrary-precision complex arithmetic for Go.
//
// It wraps the GNU MPC/MPFR/GMP libraries via cgo and exposes a Go-friendly API
// with parsing/formatting from/to strings, configurable precision, and a set of
// complex operations (algebraic + elementary + transcendental).
//
// Build requirements:
//   - libmpc, libmpfr, libgmp (headers + libs)
//     Debian/Ubuntu: sudo apt-get install -y libmpc-dev libmpfr-dev libgmp-dev build-essential
//     macOS/Homebrew: brew install mpc mpfr gmp
//
// Minimal usage:
//
//	z := apcomplex.MustParse("3.1415926535+2.718281828i", 256)
//	w := apcomplex.New(256).Exp(z)
//	fmt.Println(w.StringFixed(50))
//
// SPDX-License-Identifier: MIT
package apcomplex

/*
#cgo CFLAGS: -O2
#cgo LDFLAGS: -lmpc -lmpfr -lgmp
#include <stdlib.h>
#include <string.h>
#include <mpc.h>
#include <mpfr.h>

// Helpers to format MPFR/MPC values to C strings we can return to Go.
static char* apc_mpfr_to_str_fixed(mpfr_srcptr x, int digits) {
    if (digits < 0) digits = 0;
    int n = mpfr_snprintf(NULL, 0, "%.*Rf", digits, x);
    if (n < 0) return NULL;
    char *buf = (char*)malloc((size_t)n + 1);
    if (!buf) return NULL;
    if (mpfr_snprintf(buf, (size_t)n + 1, "%.*Rf", digits, x) < 0) {
        free(buf);
        return NULL;
    }
    return buf;
}

static char* apc_mpfr_to_str_sci(mpfr_srcptr x, int digits) {
    if (digits < 1) digits = 1;
    int n = mpfr_snprintf(NULL, 0, "%.*Re", digits, x);
    if (n < 0) return NULL;
    char *buf = (char*)malloc((size_t)n + 1);
    if (!buf) return NULL;
    if (mpfr_snprintf(buf, (size_t)n + 1, "%.*Re", digits, x) < 0) {
        free(buf);
        return NULL;
    }
    return buf;
}

static char* apc_mpc_to_a_plus_bi(mpc_srcptr z, int digits, int scientific) {
    mpfr_srcptr re = mpc_realref(z);
    mpfr_srcptr im = mpc_imagref(z);
    char *rs = scientific ? apc_mpfr_to_str_sci(re, digits) : apc_mpfr_to_str_fixed(re, digits);
    char *is = scientific ? apc_mpfr_to_str_sci(im, digits) : apc_mpfr_to_str_fixed(im, digits);
    if (!rs || !is) { if (rs) free(rs); if (is) free(is); return NULL; }
    int neg = (is[0] == '-') ? 1 : 0;
    size_t rn = strlen(rs);
    size_t in = strlen(is);
    size_t total = rn + 1 + (neg ? (in - 1) : in) + 1 + 1; // re + sign + im + 'i' + NUL
    char *out = (char*)malloc(total);
    if (!out) { free(rs); free(is); return NULL; }
    char *p = out;
    memcpy(p, rs, rn); p += rn;
    *p++ = neg ? '-' : '+';
    if (neg) { memcpy(p, is + 1, in - 1); p += in - 1; }
    else { memcpy(p, is, in); p += in; }
    *p++ = 'i';
    *p = '\0';
    free(rs); free(is);
    return out;
}

// Helpers so Go code doesn't reference MPC macros directly (cgo can't see macros).
static char* apc_mpc_real_fixed(mpc_srcptr z, int digits) {
    return apc_mpfr_to_str_fixed(mpc_realref(z), digits);
}
static char* apc_mpc_imag_fixed(mpc_srcptr z, int digits) {
    return apc_mpfr_to_str_fixed(mpc_imagref(z), digits);
}
*/
import "C"

import (
	"errors"
	"fmt"
	"runtime"
	"strings"
	"unsafe"
)

// default rounding mode (nearest, nearest)
var defaultRnd = C.mpc_rnd_t(C.MPC_RNDNN)

// Complex is an arbitrary-precision complex backed by GNU MPC/MPFR.
// Use New/Parse; zero value is not usable.
type Complex struct {
	z    C.mpc_t
	prec uint
	init bool
}

// New allocates a value with the given precision in bits (like MPFR/MPC). If bits==0, 53 is used.
func New(bits uint) *Complex {
	if bits == 0 {
		bits = 53
	}
	c := &Complex{prec: bits}
	C.mpc_init2(&c.z[0], C.mpfr_prec_t(bits))
	c.init = true
	runtime.SetFinalizer(c, func(cc *Complex) {
		if cc.init {
			C.mpc_clear(&cc.z[0])
			cc.init = false
		}
	})
	return c
}

// Close frees C resources.
func (c *Complex) Close() {
	if c != nil && c.init {
		C.mpc_clear(&c.z[0])
		c.init = false
	}
}

// Prec returns precision in bits.
func (c *Complex) Prec() uint { return c.prec }

// SetPrec changes precision (rounding value to the new precision).
func (c *Complex) SetPrec(bits uint) *Complex {
	if !c.init {
		panic("apcomplex: not initialized")
	}
	if bits == 0 {
		bits = 53
	}
	if bits == c.prec {
		return c
	}
	C.mpc_set_prec(&c.z[0], C.mpfr_prec_t(bits))
	c.prec = bits
	return c
}

// Clone returns a deep copy.
func (c *Complex) Clone() *Complex {
	out := New(c.prec)
	C.mpc_set(&out.z[0], &c.z[0], defaultRnd)
	return out
}

// Parse parses a complex literal at given precision. Accepts:
//
//	"a+bi", "a-bi", "i", "-i", plain real "a", or MPC form "(a b)" / "(a, b)".
func Parse(s string, prec uint) (*Complex, error) {
	z := New(prec)
	if err := z.SetString(s); err != nil {
		z.Close()
		return nil, err
	}
	return z, nil
}

// MustParse panics on error.
func MustParse(s string, prec uint) *Complex {
	z, err := Parse(s, prec)
	if err != nil {
		panic(err)
	}
	return z
}

// SetString sets c from a complex string (see Parse).
func (c *Complex) SetString(s string) error {
	if !c.init {
		return errors.New("apcomplex: not initialized")
	}
	re, im, ok := normalizeToPair(s)
	if !ok {
		return fmt.Errorf("apcomplex: invalid complex literal %q", s)
	}
	return c.SetBase(re, im, 10)
}

// SetBase sets c = re + i*im, parsing both parts using the given base (<=0 defaults to 10).
func (c *Complex) SetBase(re, im string, base int) error {
	if !c.init {
		return errors.New("apcomplex: not initialized")
	}
	var r, i C.mpfr_t
	C.mpfr_init2(&r[0], C.mpfr_prec_t(c.prec))
	C.mpfr_init2(&i[0], C.mpfr_prec_t(c.prec))
	defer C.mpfr_clear(&r[0])
	defer C.mpfr_clear(&i[0])

	cr := C.CString(strings.TrimSpace(re))
	ci := C.CString(strings.TrimSpace(im))
	defer C.free(unsafe.Pointer(cr))
	defer C.free(unsafe.Pointer(ci))

	b := C.int(base)
	if base <= 0 {
		b = 10
	}
	if C.mpfr_set_str(&r[0], cr, b, C.MPFR_RNDN) != 0 {
		return fmt.Errorf("apcomplex: invalid real part %q", re)
	}
	if C.mpfr_set_str(&i[0], ci, b, C.MPFR_RNDN) != 0 {
		return fmt.Errorf("apcomplex: invalid imaginary part %q", im)
	}
	C.mpc_set_fr_fr(&c.z[0], &r[0], &i[0], defaultRnd)
	return nil
}

// SetParts sets c = re + i*im using base-10 strings.
func (c *Complex) SetParts(re, im string) error { return c.SetBase(re, im, 10) }

// normalizeToPair converts common forms into separate real/imag strings.
func normalizeToPair(in string) (string, string, bool) {
	s := strings.TrimSpace(in)
	if s == "" {
		return "0", "0", true
	}
	if strings.HasPrefix(s, "(") && strings.HasSuffix(s, ")") {
		mid := strings.TrimSpace(s[1 : len(s)-1])
		mid = strings.ReplaceAll(mid, ",", " ")
		f := strings.Fields(mid)
		if len(f) == 1 {
			return f[0], "0", true
		}
		if len(f) >= 2 {
			return f[0], f[1], true
		}
		return "", "", false
	}
	s = strings.ReplaceAll(s, "I", "i")
	if s == "i" || s == "+i" {
		return "0", "1", true
	}
	if s == "-i" {
		return "0", "-1", true
	}
	if strings.HasSuffix(s, "i") {
		core := strings.TrimSpace(s[:len(s)-1])
		idx := lastSignNotInExponent(core)
		if idx > 0 {
			re := strings.TrimSpace(core[:idx])
			im := strings.TrimSpace(core[idx:])
			if im == "+" || im == "-" {
				return re, "0", true
			}
			return re, im, true
		}
		return "0", core, true
	}
	return s, "0", true
}

// lastSignNotInExponent finds last '+'/'-' not part of an exponent and not at position 0.
func lastSignNotInExponent(s string) int {
	for i := len(s) - 1; i > 0; i-- {
		if s[i] == '+' || s[i] == '-' {
			if s[i-1] != 'e' && s[i-1] != 'E' {
				return i
			}
		}
	}
	return -1
}

// Formatting
func (c *Complex) StringFixed(digits int) string {
	if digits < 0 {
		digits = 0
	}
	if !c.init {
		return "(invalid)"
	}
	p := C.apc_mpc_to_a_plus_bi(&c.z[0], C.int(digits), C.int(0))
	if p == nil {
		return "<oom>"
	}
	defer C.free(unsafe.Pointer(p))
	return C.GoString(p)
}

func (c *Complex) StringScientific(digits int) string {
	if digits < 1 {
		digits = 1
	}
	if !c.init {
		return "(invalid)"
	}
	p := C.apc_mpc_to_a_plus_bi(&c.z[0], C.int(digits), C.int(1))
	if p == nil {
		return "<oom>"
	}
	defer C.free(unsafe.Pointer(p))
	return C.GoString(p)
}

func (c *Complex) RealStringFixed(digits int) string {
	if digits < 0 {
		digits = 0
	}
	if !c.init {
		return "(invalid)"
	}
	p := C.apc_mpc_real_fixed(&c.z[0], C.int(digits))
	if p == nil {
		return "<oom>"
	}
	defer C.free(unsafe.Pointer(p))
	return C.GoString(p)
}

func (c *Complex) ImagStringFixed(digits int) string {
	if digits < 0 {
		digits = 0
	}
	if !c.init {
		return "(invalid)"
	}
	p := C.apc_mpc_imag_fixed(&c.z[0], C.int(digits))
	if p == nil {
		return "<oom>"
	}
	defer C.free(unsafe.Pointer(p))
	return C.GoString(p)
}

// Algebraic ops (mutating; return receiver for chaining)
func (c *Complex) Set(a *Complex) *Complex { C.mpc_set(&c.z[0], &a.z[0], defaultRnd); return c }
func (c *Complex) Add(a, b *Complex) *Complex {
	C.mpc_add(&c.z[0], &a.z[0], &b.z[0], defaultRnd)
	return c
}
func (c *Complex) Sub(a, b *Complex) *Complex {
	C.mpc_sub(&c.z[0], &a.z[0], &b.z[0], defaultRnd)
	return c
}
func (c *Complex) Mul(a, b *Complex) *Complex {
	C.mpc_mul(&c.z[0], &a.z[0], &b.z[0], defaultRnd)
	return c
}
func (c *Complex) Div(a, b *Complex) *Complex {
	C.mpc_div(&c.z[0], &a.z[0], &b.z[0], defaultRnd)
	return c
}
func (c *Complex) Neg(a *Complex) *Complex  { C.mpc_neg(&c.z[0], &a.z[0], defaultRnd); return c }
func (c *Complex) Conj(a *Complex) *Complex { C.mpc_conj(&c.z[0], &a.z[0], defaultRnd); return c }
func (c *Complex) Inv(a *Complex) *Complex {
	// c = 1 / a
	C.mpc_set_ui_ui(&c.z[0], 1, 0, defaultRnd)
	C.mpc_div(&c.z[0], &c.z[0], &a.z[0], defaultRnd)
	return c
}

// Elementary/transcendental
func (c *Complex) Sqrt(a *Complex) *Complex { C.mpc_sqrt(&c.z[0], &a.z[0], defaultRnd); return c }
func (c *Complex) Exp(a *Complex) *Complex  { C.mpc_exp(&c.z[0], &a.z[0], defaultRnd); return c }
func (c *Complex) Log(a *Complex) *Complex  { C.mpc_log(&c.z[0], &a.z[0], defaultRnd); return c }
func (c *Complex) Pow(a, b *Complex) *Complex {
	C.mpc_pow(&c.z[0], &a.z[0], &b.z[0], defaultRnd)
	return c
}

func (c *Complex) Sin(a *Complex) *Complex  { C.mpc_sin(&c.z[0], &a.z[0], defaultRnd); return c }
func (c *Complex) Cos(a *Complex) *Complex  { C.mpc_cos(&c.z[0], &a.z[0], defaultRnd); return c }
func (c *Complex) Tan(a *Complex) *Complex  { C.mpc_tan(&c.z[0], &a.z[0], defaultRnd); return c }
func (c *Complex) Asin(a *Complex) *Complex { C.mpc_asin(&c.z[0], &a.z[0], defaultRnd); return c }
func (c *Complex) Acos(a *Complex) *Complex { C.mpc_acos(&c.z[0], &a.z[0], defaultRnd); return c }
func (c *Complex) Atan(a *Complex) *Complex { C.mpc_atan(&c.z[0], &a.z[0], defaultRnd); return c }

func (c *Complex) Sinh(a *Complex) *Complex  { C.mpc_sinh(&c.z[0], &a.z[0], defaultRnd); return c }
func (c *Complex) Cosh(a *Complex) *Complex  { C.mpc_cosh(&c.z[0], &a.z[0], defaultRnd); return c }
func (c *Complex) Tanh(a *Complex) *Complex  { C.mpc_tanh(&c.z[0], &a.z[0], defaultRnd); return c }
func (c *Complex) Asinh(a *Complex) *Complex { C.mpc_asinh(&c.z[0], &a.z[0], defaultRnd); return c }
func (c *Complex) Acosh(a *Complex) *Complex { C.mpc_acosh(&c.z[0], &a.z[0], defaultRnd); return c }
func (c *Complex) Atanh(a *Complex) *Complex { C.mpc_atanh(&c.z[0], &a.z[0], defaultRnd); return c }

// Magnitude/argument as strings (computed with MPFR real temporaries)
func (c *Complex) AbsStringFixed(a *Complex, digits int) string {
	if digits < 0 {
		digits = 0
	}
	var r C.mpfr_t
	C.mpfr_init2(&r[0], C.mpfr_prec_t(c.prec))
	defer C.mpfr_clear(&r[0])
	C.mpc_abs(&r[0], &a.z[0], C.MPFR_RNDN)
	p := C.apc_mpfr_to_str_fixed(&r[0], C.int(digits))
	if p == nil {
		return "<oom>"
	}
	defer C.free(unsafe.Pointer(p))
	return C.GoString(p)
}

func (c *Complex) ArgStringScientific(a *Complex, digits int) string {
	if digits < 1 {
		digits = 1
	}
	var r C.mpfr_t
	C.mpfr_init2(&r[0], C.mpfr_prec_t(c.prec))
	defer C.mpfr_clear(&r[0])
	C.mpc_arg(&r[0], &a.z[0], C.MPFR_RNDN)
	p := C.apc_mpfr_to_str_sci(&r[0], C.int(digits))
	if p == nil {
		return "<oom>"
	}
	defer C.free(unsafe.Pointer(p))
	return C.GoString(p)
}

// Non-mutating convenience wrappers
func Add(a, b *Complex) *Complex { return New(a.prec).Add(a, b) }
func Sub(a, b *Complex) *Complex { return New(a.prec).Sub(a, b) }
func Mul(a, b *Complex) *Complex { return New(a.prec).Mul(a, b) }
func Div(a, b *Complex) *Complex { return New(a.prec).Div(a, b) }
func Neg(a *Complex) *Complex    { return New(a.prec).Neg(a) }
func Conj(a *Complex) *Complex   { return New(a.prec).Conj(a) }
func Inv(a *Complex) *Complex    { return New(a.prec).Inv(a) }
func Sqrt(a *Complex) *Complex   { return New(a.prec).Sqrt(a) }
func Exp(a *Complex) *Complex    { return New(a.prec).Exp(a) }
func Log(a *Complex) *Complex    { return New(a.prec).Log(a) }
func Pow(a, b *Complex) *Complex {
	p := a.prec
	if b.prec > p {
		p = b.prec
	}
	return New(p).Pow(a, b)
}
func Sin(a *Complex) *Complex   { return New(a.prec).Sin(a) }
func Cos(a *Complex) *Complex   { return New(a.prec).Cos(a) }
func Tan(a *Complex) *Complex   { return New(a.prec).Tan(a) }
func Asin(a *Complex) *Complex  { return New(a.prec).Asin(a) }
func Acos(a *Complex) *Complex  { return New(a.prec).Acos(a) }
func Atan(a *Complex) *Complex  { return New(a.prec).Atan(a) }
func Sinh(a *Complex) *Complex  { return New(a.prec).Sinh(a) }
func Cosh(a *Complex) *Complex  { return New(a.prec).Cosh(a) }
func Tanh(a *Complex) *Complex  { return New(a.prec).Tanh(a) }
func Asinh(a *Complex) *Complex { return New(a.prec).Asinh(a) }
func Acosh(a *Complex) *Complex { return New(a.prec).Acosh(a) }
func Atanh(a *Complex) *Complex { return New(a.prec).Atanh(a) }
