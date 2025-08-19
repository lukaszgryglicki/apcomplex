package apcomplex

import (
	"math"
	"math/big"
	"strconv"
	"strings"
	"testing"
)

// helper: parse with test precision
func tp(s string) *Complex { return MustParse(s, 128) }

// helper: parse decimal string (from RealStringFixed/ImagStringFixed) to float64
func f64(s string) float64 {
	// strings may be like "+0.0000" — trim leading '+'
	s = strings.TrimSpace(s)
	if len(s) > 0 && s[0] == '+' {
		s = s[1:]
	}
	v, _ := strconv.ParseFloat(s, 64)
	return v
}

// helper: |a-b| <= tol (component-wise on re & im)
func equalApprox(a, b *Complex, tol float64) bool {
	diff := Sub(a, b)
	re := f64(diff.RealStringFixed(40))
	im := f64(diff.ImagStringFixed(40))
	return math.Abs(re) <= tol && math.Abs(im) <= tol
}

func TestParseFormatRoundTrip(t *testing.T) {
	tests := []string{
		"0",
		"1",
		"-1",
		"i",
		"-i",
		"3.1415926535+2.718281828i",
		"3.1415926535-2.718281828i",
		"(2.5  -4.75)",
		"(2.5, -4.75)",
	}
	for _, s := range tests {
		z, err := Parse(s, 128)
		if err != nil {
			t.Fatalf("Parse %q failed: %v", s, err)
		}
		_ = z.StringFixed(30)      // ensure formatting works
		_ = z.StringScientific(20) // ensure sci formatting works
	}
}

func TestBasicAlgebra(t *testing.T) {
	z := tp("3.25-1.75i")
	negz := Neg(z)
	sum := Add(z, negz)
	if !equalApprox(sum, tp("0"), 1e-30) {
		t.Fatalf("z + (-z) != 0, got %s", sum.StringFixed(20))
	}

	one := tp("1")
	invz := Inv(z)
	prod := Mul(z, invz)
	if !equalApprox(prod, one, 1e-28) { // slightly looser because of division
		t.Fatalf("z * inv(z) != 1, got %s", prod.StringFixed(20))
	}

	conjz := Conj(z)
	conjConj := Conj(conjz)
	if !equalApprox(conjConj, z, 1e-30) {
		t.Fatalf("conj(conj(z)) != z, got %s vs %s", conjConj.StringFixed(20), z.StringFixed(20))
	}
}

func TestAddSubMulDiv(t *testing.T) {
	a := tp("1.5+0.75i")
	b := tp("-2.25+0.5i")

	wantAdd := tp("-0.75+1.25i")
	gotAdd := Add(a, b)
	if !equalApprox(gotAdd, wantAdd, 1e-30) {
		t.Fatalf("Add mismatch: got %s, want %s", gotAdd.StringFixed(20), wantAdd.StringFixed(20))
	}

	wantSub := tp("3.75+0.25i")
	gotSub := Sub(a, b)
	if !equalApprox(gotSub, wantSub, 1e-30) {
		t.Fatalf("Sub mismatch: got %s, want %s", gotSub.StringFixed(20), wantSub.StringFixed(20))
	}

	wantMul := tp("-3.75-0.9375i")
	gotMul := Mul(a, b)
	if !equalApprox(gotMul, wantMul, 1e-30) {
		t.Fatalf("Mul mismatch: got %s, want %s", gotMul.StringFixed(20), wantMul.StringFixed(20))
	}

	// Division check via inverse: a/b ≈ a*inv(b)
	gotDiv := Div(a, b)
	gotAlt := Mul(a, Inv(b))
	if !equalApprox(gotDiv, gotAlt, 1e-28) {
		t.Fatalf("Div mismatch a/b vs a*inv(b): %s vs %s", gotDiv.StringFixed(20), gotAlt.StringFixed(20))
	}
}

func TestExpLog(t *testing.T) {
	// avoid branch cut issues: pick a generic complex not on negative real axis
	z := tp("0.75+0.5i")
	el := Exp(Log(z)) // exp(log(z)) ~= z
	if !equalApprox(el, z, 1e-28) {
		t.Fatalf("exp(log(z)) != z, got %s vs %s", el.StringFixed(20), z.StringFixed(20))
	}
}

func TestTrigIdentityReal(t *testing.T) {
	// For real x, sin^2(x)+cos^2(x)=1. Use a real input.
	x := tp("0.5") // purely real
	s := Sin(x)
	c := Cos(x)
	s2 := Mul(s, s)
	c2 := Mul(c, c)
	sum := Add(s2, c2)
	if !equalApprox(sum, tp("1"), 1e-28) {
		t.Fatalf("sin^2+cos^2 != 1 for real x, got %s", sum.StringFixed(30))
	}
}

// --- High-precision tests for exp/log and very large integers ---

// bigPow2String returns the exact decimal string of 2^n using math/big.
func bigPow2String(n uint) string {
	b := new(big.Int).Lsh(big.NewInt(1), n)
	return b.String()
}

// trimPlusZero removes a leading '+' and normalizes "-0" to "0".
func trimPlusZero(s string) string {
	s = strings.TrimSpace(s)
	if strings.HasPrefix(s, "+") {
		s = s[1:]
	}
	if s == "-0" {
		return "0"
	}
	return s
}

func TestFormatPow2_1024_AllDigits(t *testing.T) {
	want := bigPow2String(1024)
	z, err := Parse(want, 2048)
	if err != nil {
		t.Fatalf("Parse(2^1024) failed: %v", err)
	}
	got := z.RealStringFixed(0)
	if trimPlusZero(got) != want {
		t.Fatalf("format mismatch for 2^1024: got %q", got)
	}
	if trimPlusZero(z.ImagStringFixed(0)) != "0" {
		t.Fatalf("imag part not zero: %q", z.ImagStringFixed(0))
	}
}

func TestExpLog_Pow2_1024_AllDigits(t *testing.T) {
	prec := uint(4096)
	want := bigPow2String(1024)
	two := MustParse("2", prec)
	ln2 := New(prec).Log(two)
	k := MustParse("1024", prec)
	tval := New(prec).Mul(ln2, k)
	pow := New(prec).Exp(tval) // exp(ln(2)*1024) = 2^1024
	got := pow.RealStringFixed(0)
	if trimPlusZero(got) != want {
		t.Fatalf("exp(log(2)*1024) mismatch: got %q", got)
	}
	if trimPlusZero(pow.ImagStringFixed(0)) != "0" {
		t.Fatalf("imag part of exp(log(2)*1024) not zero: %q", pow.ImagStringFixed(0))
	}
}

func TestExpLog_RoundTrip_VeryLargeComplex(t *testing.T) {
	prec := uint(1024)
	w := MustParse("1e100+1e100i", prec)
	lw := New(prec).Log(w)
	back := New(prec).Exp(lw)
	if !equalApprox(back, w, 1e-25) {
		t.Fatalf("exp(log(w)) != w for very large w, got %s vs %s",
			back.StringScientific(30), w.StringScientific(30))
	}
}

func TestLog1AndExp0Exact(t *testing.T) {
	prec := uint(256)
	one := MustParse("1", prec)
	zero := MustParse("0", prec)
	ln1 := New(prec).Log(one)
	if trimPlusZero(ln1.RealStringFixed(0)) != "0" || trimPlusZero(ln1.ImagStringFixed(0)) != "0" {
		t.Fatalf("log(1) != 0, got %s", ln1.StringFixed(0))
	}
	e0 := New(prec).Exp(zero)
	if trimPlusZero(e0.RealStringFixed(0)) != "1" || trimPlusZero(e0.ImagStringFixed(0)) != "0" {
		t.Fatalf("exp(0) != 1, got %s", e0.StringFixed(0))
	}
}
