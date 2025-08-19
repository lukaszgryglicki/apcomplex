package main

// Tetration evaluator using apcomplex (arbitrary-precision complex).
//
// Usage (exactly 3 positional args):
//   tetrate <base> <height> <precision_bits>
//
// Examples:
//   tetrate 2+1e-100i 2048 8192
//   tetrate "(0.5 0.5)" 1.5 2048
//
// Definition used:
//   Let f(z) = b^z (principal branch). We compute continuous iteration f^{∘t}(1),
//   so T_b(t) = f^{∘t}(1). For integer t, this matches the usual power tower.
//
// Method (primary): Koenigs/Schröder linearization at an attracting fixed point z*:
//   1) Find z* solving z* = b^{z*}; λ = f'(z*) = ln(b)*z*.
//   2) Koenigs map φ(z) = lim_{n→∞} λ^{-n}(f^{∘n}(z) - z*).
//   3) f^{∘t}(z) = φ^{-1}(λ^{t} φ(z)). We compute φ via iteration, and obtain φ^{-1} by
//      solving for w: φ(w) = y using Newton and an iterative approximation for φ.
//
// Fallback:
//   If |λ| ≥ 1 or z* cannot be found, we only support integer heights by direct
//   right-associated exponentiation (the classic power tower). Non-integer heights
//   in non-attracting regimes are not implemented (they require deeper Abel/Kneser
//   constructions beyond this demo).
//
// Notes:
//   * Principal complex log/exp are used (branch cuts apply).
//   * Numerical parameters (iteration counts/tolerances) are chosen from precision.
//   * This tool prints the result in scientific form with ~digits from precision,
//     plus a quick functional check value b^T(h).
//
// SPDX-License-Identifier: MIT

import (
	"errors"
	"flag"
	"fmt"
	"math"
	"os"
	"strconv"
	"strings"

	ap "github.com/lukaszgryglicki/apcomplex"
)

func main() {
	flag.CommandLine.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage: %s <base> <height> <precision_bits>\n", os.Args[0])
		fmt.Fprintf(os.Stderr, "Example: %s 2+1e-100i 2048 8192\n", os.Args[0])
	}
	flag.Parse()
	if flag.NArg() != 3 {
		flag.Usage()
		os.Exit(2)
	}

	baseStr := flag.Arg(0)
	heightStr := flag.Arg(1)
	bits64, err := strconv.ParseUint(flag.Arg(2), 10, 32)
	if err != nil || bits64 == 0 {
		fmt.Fprintln(os.Stderr, "invalid precision bits; need positive integer")
		os.Exit(2)
	}
	prec := uint(bits64)

	b, err := ap.Parse(baseStr, prec)
	if err != nil {
		fmt.Fprintln(os.Stderr, "parse base:", err)
		os.Exit(2)
	}
	h, err := ap.Parse(heightStr, prec)
	if err != nil {
		fmt.Fprintln(os.Stderr, "parse height:", err)
		os.Exit(2)
	}

	res, method, err := tetrate(b, h, prec)
	if err != nil {
		fmt.Fprintln(os.Stderr, "tetrate:", err)
		os.Exit(1)
	}

	digits := int(math.Max(1, math.Floor(float64(prec)*math.Log10(2))-5))
	fmt.Printf("method: %s\n", method)
	fmt.Printf("T_b(h) with b=%s, h=%s\n", b.StringScientific(digits), h.StringScientific(digits))
	fmt.Printf("precision=%d bits, printing ~%d significant digits\n", prec, digits)
	fmt.Printf("result: %s\n", res.StringScientific(digits))
	fmt.Printf("Re: %s\n", res.RealStringFixed(0))
	fmt.Printf("Im: %s\n", res.ImagStringFixed(0))
	// quick functional check: b^{T(h)} should equal T(h+1)
	bpow := ap.New(prec)
	bpow.Exp(ap.New(prec).Mul(ap.New(prec).Log(b), res)) // e^{log(b)*res} = b^res
	fmt.Printf("b^(T(h)) (sanity): %s\n", bpow.StringScientific(digits))
}

// tetrate computes T_b(h) = f^{∘h}(1) where f(z)=b^z.
// Returns result, method description, error.
func tetrate(b, h *ap.Complex, prec uint) (*ap.Complex, string, error) {
	// Special bases
	if isApproximatelyOne(b, prec) {
		// f(z)=1^z = 1; then f^{∘t}(1)=1 for all t.
		return ap.MustParse("1", prec), "constant base=1", nil
	}

	f := func(z *ap.Complex) *ap.Complex {
		return ap.New(prec).Exp(ap.New(prec).Mul(ap.New(prec).Log(b), z)) // b^z
	}

	// Try Schröder/Koenigs in basin of attraction.
	zstar, lam, ok := findAttractingFixedPoint(b, f, prec)
	if ok {
		res, err := tetrateSchroeder(b, f, zstar, lam, h, prec)
		if err == nil {
			return res, "Schröder (Koenigs) fractional iteration", nil
		}
		// fall through to possible integer fallback
	}

	// Fallback: if height is a non-negative integer, compute classic tower.
	if n, isInt := tryIntegerHeight(h, prec); isInt && n >= 0 {
		res := powerTower(b, n, prec)
		return res, "integer tower fallback", nil
	}

	return nil, "", errors.New("non-attracting regime or fixed point not found; non-integer heights require advanced Abel/Kneser methods not implemented here")
}

// tetrateSchroeder computes f^{∘h}(1) via Koenigs map φ and Newton inversion.
func tetrateSchroeder(b *ap.Complex, f func(*ap.Complex) *ap.Complex, zstar, lam, h *ap.Complex, prec uint) (*ap.Complex, error) {
	// Choose K such that |λ|^K is tiny relative to target precision.
	lamAbs := absFloat(lam, prec)
	if lamAbs <= 0 {
		return nil, errors.New("invalid derivative magnitude at fixed point")
	}

	// desired decimal digits
	digs := float64(prec)*math.Log10(2)
	// target error for φ ~ 10^{-digs/2}
	var K int
	if lamAbs < 1 {
		K = int(math.Ceil((digs/2.0)/(-math.Log10(lamAbs))))
		if K < 8 {
			K = 8
		}
		if K > 2000 {
			K = 2000 // hard cap
		}
	} else {
		// shouldn't happen if attracting
		return nil, errors.New("|λ| >= 1 (not attracting)")
	}

	// φ(1)
	phi1 := koenigsPhi(f, zstar, lam, ap.MustParse("1", prec), K, prec)
	// λ^h
	lnLam := ap.New(prec).Log(lam)
	lamPowH := ap.New(prec).Exp(ap.New(prec).Mul(lnLam, h))
	// target y = λ^h * φ(1)
	y := ap.New(prec).Mul(lamPowH, phi1)

  // Solve φ(w) = y for w via Newton on Φ_K(w) = λ^{-K}(f^{∘K}(w) - z*)
  w0 := ap.New(prec).Add(zstar, y) // first-order inverse near z*
  w := newtonSolvePhiEqWithLnB(f, zstar, lam, y, w0, K, prec)
	return w, nil
}

// findAttractingFixedPoint tries to locate z* = b^{z*} with |λ|<1 where λ = ln(b)*z*.
func findAttractingFixedPoint(b *ap.Complex, f func(*ap.Complex) *ap.Complex, prec uint) (zstar, lam *ap.Complex, ok bool) {
	// iterate from 1: u_{n+1} = f(u_n)
	u := ap.MustParse("1", prec)
	var last *ap.Complex
	for i := 0; i < 2000; i++ {
		last = u
		u = f(u)
		if diffSmall(u, last, prec, 20) { // converged
			zstar = u
			lnb := ap.New(prec).Log(b)
			lam = ap.New(prec).Mul(lnb, zstar)
			if absFloat(lam, prec) < 1 {
				return zstar, lam, true
			}
			break
		}
		if absFloat(u, prec) > 1e12 { // likely diverging
			break
		}
	}
	// Try Newton on g(z)=z - f(z)
	g := func(z *ap.Complex) *ap.Complex { return ap.New(prec).Sub(z, f(z)) }
	gp := func(z *ap.Complex) *ap.Complex {
		// g'(z) = 1 - f'(z) = 1 - ln(b)*b^z
		one := ap.MustParse("1", prec)
		lnb := ap.New(prec).Log(b)
		fz := f(z)
		fp := ap.New(prec).Mul(lnb, fz)
		return ap.New(prec).Sub(one, fp)
	}
	z0 := ap.MustParse("1", prec)
	z := z0
	for i := 0; i < 100; i++ {
		gz := g(z)
		gprime := gp(z)
		step := ap.New(prec).Div(gz, gprime)
		z = ap.New(prec).Sub(z, step)
		if diffSmall(gz, ap.MustParse("0", prec), prec, 30) {
			break
		}
	}
	zstar = z
	lnb := ap.New(prec).Log(b)
	lam = ap.New(prec).Mul(lnb, zstar)
	if absFloat(lam, prec) < 1 {
		return zstar, lam, true
	}
	return nil, nil, false
}

// koenigsPhi approximates φ(z) ≈ λ^{-K}(f^{∘K}(z) - z*).
func koenigsPhi(f func(*ap.Complex) *ap.Complex, zstar, lam, z *ap.Complex, K int, prec uint) *ap.Complex {
	u := z
	for i := 0; i < K; i++ { u = f(u) }
	num := ap.New(prec).Sub(u, zstar)
	// lamPow := ap.New(prec).Exp(ap.New(prec).Mul(ap.New(prec).Log(lam), ap.MustParse(strconvI(-K), prec)))
	// but more stable to compute λ^{-K} as (λ^K)^{-1}:
	lamPowPos := ap.New(prec).Exp(ap.New(prec).Mul(ap.New(prec).Log(lam), ap.MustParse(strconvI(K), prec)))
	lamInvPow := ap.New(prec).Inv(lamPowPos)
	return ap.New(prec).Mul(lamInvPow, num)
}

// We wrap the iterate with explicit ln(b) to compute derivatives.
func newtonSolvePhiEqWithLnB(f func(*ap.Complex) *ap.Complex, zstar, lam, y, w0 *ap.Complex, K int, prec uint) *ap.Complex {
	// Recover ln(b) from f by probing at 1: f(1)=b^1=b so ln(b)=log(f(1)) - 0
	b := f(ap.MustParse("1", prec))
	lnb := ap.New(prec).Log(b)
	w := w0
	one := ap.MustParse("1", prec)
	lamPowPos := ap.New(prec).Exp(ap.New(prec).Mul(ap.New(prec).Log(lam), ap.MustParse(strconvI(K), prec)))
	lamInvPow := ap.New(prec).Inv(lamPowPos)
	for it := 0; it < 80; it++ {
		// forward K iterations and derivative
		u := w
		der := one
		for k := 0; k < K; k++ {
			v := ap.New(prec).Exp(ap.New(prec).Mul(lnb, u)) // f(u)
			fp := ap.New(prec).Mul(lnb, v)                 // f'(u)
			der = ap.New(prec).Mul(der, fp)
			u = v
		}
		phiApprox := ap.New(prec).Mul(lamInvPow, ap.New(prec).Sub(u, zstar))
		resid := ap.New(prec).Sub(phiApprox, y)
		if diffSmall(resid, ap.MustParse("0", prec), prec, 40) { return w }
		// Newton step: w -= (phiApprox - y)/phiApprox'
		phiPrime := ap.New(prec).Mul(lamInvPow, der)
		step := ap.New(prec).Div(resid, phiPrime)
		w = ap.New(prec).Sub(w, step)
	}
	return w
}

// powerTower computes f^{∘n}(1) for integer n >= 0 (right-associated tower).
func powerTower(b *ap.Complex, n int, prec uint) *ap.Complex {
	x := ap.MustParse("1", prec)
	for i := 0; i < n; i++ {
		x = ap.New(prec).Exp(ap.New(prec).Mul(ap.New(prec).Log(b), x))
	}
	return x
}

// tryIntegerHeight returns (n, true) if h is a real integer within small tolerance.
func tryIntegerHeight(h *ap.Complex, prec uint) (int, bool) {
	// Check imag ~ 0 and real ~ integer
	im := h.ImagStringFixed(0)
	if strings.TrimSpace(im) != "0" { return 0, false }
	reStr := h.RealStringFixed(0)
	n, err := strconv.Atoi(strings.TrimSpace(strings.TrimPrefix(reStr, "+")))
	if err != nil { return 0, false }
	return n, true
}

// Helpers ------------------------------------------------------------

func diffSmall(a, b *ap.Complex, prec uint, digits int) bool {
	d := ap.New(prec).Sub(a, b)
	mag := absFloat(d, prec)
	// threshold ~ 10^{-digits}
	thr := math.Pow10(-digits) // Go 1.20+: beware Pow10 only for ints; implement via Pow(10,-digits)
	thr = math.Pow(10, float64(-digits))
	return mag < thr
}

func absFloat(a *ap.Complex, prec uint) float64 {
	tmp := ap.New(prec)
	s := tmp.AbsStringFixed(a, 30)
	val, _ := strconv.ParseFloat(strings.TrimSpace(strings.TrimPrefix(s, "+")), 64)
	return val
}

func isApproximatelyOne(b *ap.Complex, prec uint) bool {
	ones := ap.MustParse("1", prec)
	d := ap.New(prec).Sub(b, ones)
	return absFloat(d, prec) < 1e-30
}

func strconvI(i int) string { return strconv.Itoa(i) }

